# Plate buckling benchmarks for TriShellFiniteElement.jl
#
# Rectangular plate a × b × t in uniform compression along a (loading direction ξ, width direction η).
# Classical results:  σ_cr = k π² E / (12 (1-ν²)) (t/b)²
#   simply supported all edges, integer a/b : k = 4.0
#   clamped all edges, square, uniaxial      : k ≈ 10.07  (Timoshenko & Gere)
#
# The pre-buckling membrane stress is computed with the element itself (tributary edge loads,
# symmetry-line datums so the stress field is exactly uniform), then K φ = P_cr (-Kg) φ is solved
# densely on the free dofs.

using TriShellFiniteElement, Ferrite, LinearAlgebra, Test

const E = 200000.0
const ν = 0.30

# Structured triangle mesh. plane = :xy puts the plate in the global XY plane (normal Z);
# plane = :yz makes ξ = Z and η = Y (normal X), like a stud web.
function plate_grid(a, b, nx, ny; plane = :xy)
    nodes = Ferrite.Node{3,Float64}[]
    for j in 0:ny, i in 0:nx
        ξ = a * i / nx; η = b * j / ny
        push!(nodes, Ferrite.Node(plane == :xy ? Ferrite.Vec(ξ, η, 0.0) : Ferrite.Vec(0.0, η, ξ)))
    end
    id(i, j) = j * (nx + 1) + i + 1
    cells = Ferrite.Triangle[]
    for j in 0:ny-1, i in 0:nx-1
        n1, n2, n3, n4 = id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)
        push!(cells, Ferrite.Triangle((n1, n2, n3)))
        push!(cells, Ferrite.Triangle((n1, n3, n4)))
    end
    return Ferrite.Grid(cells, nodes)
end

# (ξ, η, normal) component indices in the :u field and the (ξ, η) coordinates of a node
comp(plane) = plane == :xy ? (1, 2, 3) : (3, 2, 1)
ξη(node, plane) = plane == :xy ? (node.x[1], node.x[2]) : (node.x[3], node.x[2])

"""
Return the buckling coefficient k of the first mode and the axial shortening under the unit load.
"""
function plate_buckling(a, b, t, nx, ny; Cs = TriShellFiniteElement.DEFAULT_SHEAR_RELAXATION,
                        plane = :xy, clamped = false)
    grid = plate_grid(a, b, nx, ny; plane)
    ip  = Lagrange{RefTriangle,1}()
    ip3 = TriShellFiniteElement.IP3(); ip6 = TriShellFiniteElement.IP6()
    qr1 = QuadratureRule{RefTriangle}(1); qr3 = QuadratureRule{RefTriangle}(2)
    dh = DofHandler(grid); add!(dh, :u, ip^3); add!(dh, :θ, ip^3); close!(dh)

    node_to_dofs = Dict{Int,Vector{Int}}()
    for cell in CellIterator(dh)
        cd = celldofs(cell)
        for (i, n) in enumerate(cell.nodes)
            node_to_dofs[n] = [cd[3(i-1)+1], cd[3(i-1)+2], cd[3(i-1)+3], cd[9+3(i-1)+1], cd[9+3(i-1)+2], cd[9+3(i-1)+3]]
        end
    end

    K = allocate_matrix(dh)
    K = TriShellFiniteElement.assemble_global_Ke!(K, dh, qr1, qr3, ip3, ip6, E, ν, t; Cs = Cs)

    cξ, cη, cn = comp(plane)
    tol = 1e-8
    nodes_where(f) = Set(i for (i, n) in enumerate(grid.nodes) if f(ξη(n, plane)...))
    edge0   = nodes_where((ξ, η) -> abs(ξ) < tol)
    edgeA   = nodes_where((ξ, η) -> abs(ξ - a) < tol)
    edgeB0  = nodes_where((ξ, η) -> abs(η) < tol)
    edgeBb  = nodes_where((ξ, η) -> abs(η - b) < tol)
    midline = nodes_where((ξ, η) -> abs(ξ - a / 2) < tol)   # u_ξ datum (symmetric loading)
    centerl = nodes_where((ξ, η) -> abs(η - b / 2) < tol)   # u_η datum (symmetric Poisson expansion)
    @assert !isempty(midline) && !isempty(centerl) "nx and ny must be even"
    all_edges = union(edge0, edgeA, edgeB0, edgeBb)

    # unit total compression at each loaded edge, tributary-length weighted (self-equilibrated)
    F = zeros(ndofs(dh))
    for (edge, sgn) in ((edgeA, -1.0), (edge0, +1.0)), n in edge
        (ξ, η) = ξη(grid.nodes[n], plane)
        w = (abs(η) < tol || abs(η - b) < tol) ? 0.5 : 1.0
        F[node_to_dofs[n][cξ]] += sgn * w / ny
    end

    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, all_edges, (x, t_) -> [0.0], [cn]))                  # w = 0 on all edges
    clamped && add!(ch, Dirichlet(:θ, all_edges, (x, t_) -> [0.0, 0.0, 0.0], [1, 2, 3]))
    add!(ch, Dirichlet(:u, centerl, (x, t_) -> [0.0], [cη]))
    add!(ch, Dirichlet(:u, midline, (x, t_) -> [0.0], [cξ]))
    close!(ch)

    apply!(K, F, ch)
    u = K \ F
    apply!(u, ch)
    shortening = sum(u[node_to_dofs[n][cξ]] for n in edge0) / length(edge0) -
                 sum(u[node_to_dofs[n][cξ]] for n in edgeA) / length(edgeA)

    # membrane stresses in each element's local frame
    D_stress = TriShellFiniteElement.calculate_membrane_constitutive_matrix(E, ν, t) / t
    cv = CellValues(qr1, ip3, ip3)
    σXX = zeros(getncells(grid)); σYY = zeros(getncells(grid)); τXY = zeros(getncells(grid))
    for cell in CellIterator(dh)
        xg = getcoordinates(cell)
        T = TriShellFiniteElement.calculation_rotation_matrix(xg)
        xl = TriShellFiniteElement.global_nodal_coords_to_planar_coords(xg, T)
        reinit!(cv, xl)
        dNdx = cv.fun_values.dNdx
        B = hcat([[dNdx[i][1] 0.0; 0.0 dNdx[i][2]; dNdx[i][2] dNdx[i][1]] for i in 1:3]...)
        ul = Float64[]
        for n in cell.nodes
            append!(ul, (T' * [u[node_to_dofs[n][k]] for k in 1:3])[1:2])
        end
        σ = D_stress * B * ul
        c = cellid(cell); σXX[c] = σ[1]; σYY[c] = σ[2]; τXY[c] = σ[3]
    end

    Kg = allocate_matrix(dh)
    Kg = TriShellFiniteElement.assemble_global_Kg!(Kg, dh, qr1, ip3, σXX .* t, σYY .* t, τXY .* t)

    free = setdiff(1:ndofs(dh), ch.prescribed_dofs)
    Kff  = Symmetric(Matrix(K[free, free]))
    Kgff = Symmetric(Matrix(Kg[free, free]))
    μ = eigvals(-Kgff, Kff)                    # μ = 1 / P  for  K φ = P (-Kg) φ
    Pcr = 1 / maximum(μ)
    k = Pcr / (b * t) * 12 * (1 - ν^2) * b^2 / (π^2 * E * t^2)
    return (; k, Pcr, shortening)
end

b = 92.1              # mm (362S162 web flat width)
t = 0.0346 * 25.4     # mm (33 mil), b/t ≈ 105

@testset "TriShellFiniteElement" begin

    @testset "element matrices" begin
        x = [Ferrite.Vec((0.0, 0.0)), Ferrite.Vec((10.0, 0.0)), Ferrite.Vec((0.0, 10.0))]
        ip3 = TriShellFiniteElement.IP3(); ip6 = TriShellFiniteElement.IP6()
        qr1 = QuadratureRule{RefTriangle}(1); qr3 = QuadratureRule{RefTriangle}(2)
        ke0 = TriShellFiniteElement.local_elastic_stiffness_matrix!(qr1, qr3, ip3, ip6, E, ν, 1.0, x; Cs = 0.0)
        ke  = TriShellFiniteElement.local_elastic_stiffness_matrix!(qr1, qr3, ip3, ip6, E, ν, 1.0, x)
        @test size(ke) == (18, 18)
        @test ke0 ≈ ke0' && ke ≈ ke'                    # symmetric to round-off (static condensation)
        @test all(diag(ke) .> 0)
        # relaxation reduces (never increases) stiffness and leaves the membrane block untouched
        @test all(eigvals(Symmetric(ke0 - ke)) .> -1e-8 * maximum(abs, ke0))
        m = [1, 2, 7, 8, 13, 14]
        @test ke[m, m] ≈ ke0[m, m]
        @test TriShellFiniteElement.DEFAULT_SHEAR_RELAXATION == 0.2
    end

    @testset "membrane: uniform compression" begin
        r = plate_buckling(b, b, t, 6, 6)
        @test r.shortening ≈ 1.0 * b / (E * b * t) rtol = 1e-8      # P L / (E A)
    end

    @testset "simply supported plate, k = 4 (Cs = 0.2 default)" begin
        @test plate_buckling(b, b, t, 10, 10).k ≈ 4.0 rtol = 0.02
        @test plate_buckling(b, b, t, 16, 16).k ≈ 4.0 rtol = 0.01
        @test plate_buckling(3b, b, t, 30, 10).k ≈ 4.0 rtol = 0.02   # a/b = 3, three half-waves
    end

    @testset "unrelaxed element (Cs = 0) locks on a coarse mesh" begin
        k0 = plate_buckling(b, b, t, 10, 10; Cs = 0.0).k
        @test k0 ≈ 4.4392 rtol = 1e-3                                # documented legacy value
        @test plate_buckling(b, b, t, 24, 24; Cs = 0.0).k ≈ 4.0 rtol = 0.01   # converges with refinement
        @test plate_buckling(b, b, 0.2, 10, 10; Cs = 0.0).k > 8.0    # thin plate: locking blows up
    end

    @testset "thickness independence with Cs = 0.2" begin
        # b/t = 460 … 46; the thickest case carries a genuine ~0.3% Mindlin shear reduction
        ks = [plate_buckling(b, b, tt, 10, 10).k for tt in (0.2, 0.5, t, 2.0)]
        @test all(isapprox.(ks, 4.0; rtol = 0.02))
        @test maximum(ks) - minimum(ks) < 0.10
        # the unrelaxed element on the same meshes spans k = 4.1 … 8.3
        @test plate_buckling(b, b, 0.2, 10, 10; Cs = 0.0).k - plate_buckling(b, b, 2.0, 10, 10; Cs = 0.0).k > 4.0
    end

    @testset "orientation invariance in 3D" begin
        kxy = plate_buckling(b, b, t, 8, 8; plane = :xy).k
        kyz = plate_buckling(b, b, t, 8, 8; plane = :yz).k
        @test kxy ≈ kyz rtol = 1e-8
    end

    @testset "clamped square plate, k ≈ 10.07" begin
        @test plate_buckling(b, b, t, 16, 16; clamped = true).k ≈ 10.07 rtol = 0.04
    end
end
