# Tests of the elastic–plastic, geometrically nonlinear capability (plasticity.jl, nonlinear.jl,
# solver.jl). Included by runtests.jl.

using TriShellFiniteElement, Ferrite, LinearAlgebra, StaticArrays, Random, Test
const TS = TriShellFiniteElement

# true stress – plastic strain table of Moen (2008) Fig. 7.1 (Yu 2005 coupon), ksi
const HARD = TS.IsotropicHardening([0.0, 0.0016, 0.0259, 0.07, 0.1407, 0.1988], [33.0, 34.31, 42.14, 51.35, 60.86, 65.5])
const MAT = TS.J2PlaneStress(29500.0, 0.3, HARD)

function tri_plate_grid(a, b, nx, ny)
    nodes = Ferrite.Node{3,Float64}[]
    for j in 0:ny, i in 0:nx
        push!(nodes, Ferrite.Node(Ferrite.Vec(a * i / nx, b * j / ny, 0.0)))
    end
    id(i, j) = j * (nx + 1) + i + 1
    cells = Ferrite.Triangle[]
    for j in 0:ny-1, i in 0:nx-1
        n1, n2, n3, n4 = id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)
        push!(cells, Ferrite.Triangle((n1, n2, n3))); push!(cells, Ferrite.Triangle((n1, n3, n4)))
    end
    return Ferrite.Grid(cells, nodes)
end

function tri_dofhandler(grid)
    ip = Lagrange{RefTriangle,1}()
    dh = DofHandler(grid); add!(dh, :u, ip^3); add!(dh, :θ, ip^3); close!(dh)
    return dh
end

node_dofs(dh, n, comp) = begin
    for cell in CellIterator(dh)
        for (i, nn) in enumerate(cell.nodes)
            nn == n && return celldofs(cell)[comp <= 3 ? 3(i-1)+comp : 9+3(i-1)+comp-3]
        end
    end
    error("node not found")
end

@testset "plane-stress return mapping" begin
    # hardening curve interpolation
    @test TS.yield_stress(HARD, 0.0) == (33.0, (34.31 - 33.0) / 0.0016)
    @test TS.yield_stress(HARD, 0.0259)[1] ≈ 42.14
    @test TS.yield_stress(HARD, 0.5) == (65.5, 0.0)
    @test TS.yield_stress(TS.IsotropicHardening(50.0), 0.1) == (50.0, 0.0)
    @test TS.yield_stress(TS.IsotropicHardening(50.0, 1000.0), 0.01) == (60.0, 1000.0)

    # elastic step returns the elastic tangent
    σ, C, st = TS.return_mapping(MAT, SVector(5e-4, -1.5e-4, 0.0), TS.LayerState())
    @test C == MAT.C
    @test st.ε̄p == 0.0
    @test σ ≈ MAT.C * SVector(5e-4, -1.5e-4, 0.0)

    # plastic step: stress on the hardened yield surface, consistent tangent = finite difference
    st0 = TS.LayerState(SVector(0.004, -0.002, 0.001), 0.004)
    ε0 = SVector(0.006, -0.001, 0.003)
    σ0, Cep, st1 = TS.return_mapping(MAT, ε0, st0)
    @test TS.von_mises(σ0) ≈ TS.yield_stress(HARD, st1.ε̄p)[1] rtol = 1e-9
    @test st1.ε̄p > st0.ε̄p
    @test Cep ≈ Cep' atol = 1e-9
    h = 1e-7
    for j in 1:3
        e = SVector{3}(i == j ? h : 0.0 for i in 1:3)
        col = (TS.return_mapping(MAT, ε0 + e, st0)[1] - TS.return_mapping(MAT, ε0 - e, st0)[1]) / (2h)
        @test col ≈ Cep[:, j] rtol = 1e-6
    end

    # uniaxial loading: σ = σy(ε̄p) and ε = σ/E + ε̄p along the whole curve
    st = TS.LayerState(); ε = SVector(0.0, 0.0, 0.0)
    for εx in range(0.0005, 0.05, 100)
        ε = SVector(εx, ε[2] + (εx - ε[1]) * (st.ε̄p > 0 ? -0.5 : -0.3), 0.0)
        local σ, C, s
        for it in 1:30
            σ, C, s = TS.return_mapping(MAT, ε, st)
            r = SVector(σ[2], σ[3])
            norm(r) < 1e-10 && break
            J = SMatrix{2,2}(C[2, 2], C[3, 2], C[2, 3], C[3, 3])
            Δ = J \ r
            ε = SVector(ε[1], ε[2] - Δ[1], ε[3] - Δ[2])
        end
        st = s
        @test abs(εx - (σ[1] / 29500 + st.ε̄p)) < 1e-9
        st.ε̄p > 0 && @test σ[1] ≈ TS.yield_stress(HARD, st.ε̄p)[1] rtol = 1e-8
    end
    @test st.ε̄p > 0.048

    # layers
    z, w = TS.thickness_layers(2.0, 5; rule = :simpson)
    @test z == [-1.0, -0.5, 0.0, 0.5, 1.0] && sum(w) ≈ 2.0 && sum(w .* z .^ 2) ≈ 2 / 3
    z, w = TS.thickness_layers(2.0, 5)
    @test sum(w) ≈ 2.0 && sum(w .* z .^ 2) ≈ 2 / 3
end

@testset "nonlinear element: linear limit and consistent tangent" begin
    E = 29500.0; ν = 0.3; t = 0.0346
    grid = tri_plate_grid(2.0, 1.0, 4, 2)
    dh = tri_dofhandler(grid)
    sec = TS.ShellSection(MAT, t; nlayers = 5)
    model = TS.PlasticShellModel(dh, sec)
    K = allocate_matrix(dh); f = zeros(ndofs(dh))
    TS.assemble_internal!(K, f, model, zeros(ndofs(dh)))
    ip3 = TS.IP3(); ip6 = TS.IP6(); qr1 = QuadratureRule{RefTriangle}(1); qr3 = QuadratureRule{RefTriangle}(2)
    Klin = allocate_matrix(dh); TS.assemble_global_Ke!(Klin, dh, qr1, qr3, ip3, ip6, E, ν, t)
    @test norm(f) == 0
    @test norm(K - Klin) / norm(Klin) < 1e-12      # tangent at zero = linear element (same Cs, drilling)

    # finite-difference check of the global tangent in the elastic and plastic ranges
    function fd_err(u; h = 1e-7)
        K = allocate_matrix(dh); f = zeros(ndofs(dh))
        TS.assemble_internal!(K, f, model, u)
        fp = zeros(length(u)); fm = zeros(length(u)); Kt = allocate_matrix(dh)
        err = 0.0; ref = 0.0
        for j in round.(Int, range(1, length(u), 15))
            up = copy(u); up[j] += h; TS.assemble_internal!(Kt, fp, model, up)
            um = copy(u); um[j] -= h; TS.assemble_internal!(Kt, fm, model, um)
            col = (fp - fm) / (2h)
            err = max(err, norm(K[:, j] - col)); ref = max(ref, norm(col))
        end
        return err / ref
    end
    Random.seed!(1)
    u1 = 1e-4 * randn(ndofs(dh)); u1[1:3*getnnodes(grid)] .*= 0.1
    @test fd_err(u1) < 1e-8
    u2 = 5e-2 * randn(ndofs(dh)); u2[1:3*getnnodes(grid)] .*= 0.1
    @test fd_err(u2) < 1e-6
    TS.assemble_internal!(K, f, model, u2); TS.commit!(model)
    @test maximum(r.ε̄p for r in TS.element_results(model, u2)) > 0.01   # the check was in the plastic range
end

@testset "strip in compression follows the hardening curve" begin
    t = 0.0346; a = 4.0; b = 1.0
    grid = tri_plate_grid(a, b, 8, 2)
    dh = tri_dofhandler(grid)
    model = TS.PlasticShellModel(dh, TS.ShellSection(MAT, t; nlayers = 5))
    tol = 1e-8
    left = Set(i for (i, n) in enumerate(grid.nodes) if abs(n.x[1]) < tol)
    right = Set(i for (i, n) in enumerate(grid.nodes) if abs(n.x[1] - a) < tol)
    corner = Set(i for (i, n) in enumerate(grid.nodes) if abs(n.x[1]) < tol && abs(n.x[2]) < tol)
    allnodes = Set(1:getnnodes(grid))
    δmax = 0.04 * a
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, left, (x, t_) -> [0.0], [1]))
    add!(ch, Dirichlet(:u, right, (x, t_) -> [-t_ * δmax], [1]))
    add!(ch, Dirichlet(:u, corner, (x, t_) -> [0.0], [2]))
    add!(ch, Dirichlet(:u, allnodes, (x, t_) -> [0.0], [3]))
    add!(ch, Dirichlet(:θ, allnodes, (x, t_) -> [0.0, 0.0, 0.0], [1, 2, 3]))
    close!(ch)
    rdofs = unique([node_dofs(dh, n, 1) for n in right])
    hist = TS.solve_nonlinear(model, ch; reaction_dofs = rdofs, Δt = 0.02, Δt_max = 0.05)
    @test hist.converged
    @test maximum(hist.iterations) <= 4
    for (tt, P) in zip(hist.t, hist.P)
        ε = -tt * δmax / a                 # nominal strain
        εGL = ε + ε^2 / 2                  # Green–Lagrange strain
        S = P / (b * t) / (1 + ε)          # 2nd Piola–Kirchhoff stress from the nominal reaction
        εp = -(εGL - S / 29500.0)
        εp > 1e-6 && @test -S ≈ TS.yield_stress(HARD, εp)[1] rtol = 1e-6
    end
    @test -hist.P[end] / (b * t) ≈ 42.82 rtol = 1e-3           # (1 + ε) σy(ε̄p = 0.0377)
    res = TS.element_results(model, hist.u)
    @test all(r.ε̄p ≈ res[1].ε̄p for r in res)                   # uniform state
end

@testset "elastic–plastic bending: pure bending of a strip by end rotations" begin
    # strip bent by equal and opposite end rotations, ends free to approach (no axial force), v = 0
    # and θx = 0 (plane strain across the width): uniform curvature κ = 8 κy, perfectly plastic.
    # The fully plastic plane-strain moment is M_p = (2/√3) σy t² / 4; at κ = 8 κy the exact
    # rectangular section gives M/M_p = 1 − 1/(3·64) = 0.995 and 9 Simpson layers 0.9994.
    E = 29500.0; ν = 0.3; t = 0.1; σy = 33.0
    mat = TS.J2PlaneStress(E, ν, σy)
    a = 2.0; b = 0.5
    grid = tri_plate_grid(a, b, 8, 2)
    dh = tri_dofhandler(grid)
    sec = TS.ShellSection(mat, t; nlayers = 9, rule = :simpson)
    model = TS.PlasticShellModel(dh, sec)
    tol = 1e-8
    allnodes = Set(1:getnnodes(grid))
    x0 = Set(i for (i, n) in enumerate(grid.nodes) if abs(n.x[1]) < tol)
    xa = Set(i for (i, n) in enumerate(grid.nodes) if abs(n.x[1] - a) < tol)
    κy = 2 * σy / (E * t); κmax = 8 * κy
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, allnodes, (x, t_) -> [0.0], [2]))
    add!(ch, Dirichlet(:θ, allnodes, (x, t_) -> [0.0, 0.0], [1, 3]))
    add!(ch, Dirichlet(:u, x0, (x, t_) -> [0.0, 0.0], [1, 3]))
    add!(ch, Dirichlet(:u, xa, (x, t_) -> [0.0], [3]))
    add!(ch, Dirichlet(:θ, x0, (x, t_) -> [t_ * κmax * a / 2], [2]))      # θy = −κ (x − a/2)
    add!(ch, Dirichlet(:θ, xa, (x, t_) -> [-t_ * κmax * a / 2], [2]))
    close!(ch)
    hist = TS.solve_nonlinear(model, ch; Δt = 0.05, Δt_max = 0.05)
    @test hist.converged
    res = TS.element_results(model, hist.u)
    Mx = [res[c].M[1] for c in 1:2:length(res)]          # elements whose local x axis is the global x axis
    Nx = [res[c].N[1] for c in 1:2:length(res)]
    Mp = 2 / sqrt(3) * σy * t^2 / 4
    @test all(abs.(Mx .- Mx[1]) .< 1e-4 * abs(Mx[1]))           # uniform moment
    @test -Mx[1] ≈ 0.9994 * Mp rtol = 2e-3
    @test maximum(abs, Nx) < 1e-6 * σy * t                     # no axial force
    @test maximum(r.ε̄p for r in res) > 0.005
end
