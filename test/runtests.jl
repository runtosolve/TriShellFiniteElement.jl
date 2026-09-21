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
                        plane = :xy, clamped = false, drilling = TriShellFiniteElement.DEFAULT_DRILLING)
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
    K = TriShellFiniteElement.assemble_global_Ke!(K, dh, qr1, qr3, ip3, ip6, E, ν, t; Cs, drilling)

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

# --------------------------------------------------------------------------------------------
# Saint-Venant torsion of a prismatic open section by a warping-free static twist (Moen 2008,
# Sec. 4.2.7.3.2.3, Fig. 4.41): z = 0 all nodes fixed in X, Y and one node fixed in Z; z = L all nodes
# given the rigid rotation βo about the section's mean point; T = resultant end moment, J = T L/(G βo).
# `section` is a centerline polyline (X, Y); each quad of the extruded mesh is split into two triangles.
# --------------------------------------------------------------------------------------------
function strip_section(α_deg, b; n = 8)
    a = deg2rad(α_deg); h = n ÷ 2
    X = [b / 2 * k / h for k in 0:h]; Y = zeros(h + 1)
    append!(X, [b / 2 + b / 2 * cos(a) * k / h for k in 1:h]); append!(Y, [b / 2 * sin(a) * k / h for k in 1:h])
    return X, Y
end

function lipped_c_section(; lip = 19.05, B = 76.2, D = 76.2, R = 5.05, n_flat = 4, n_corner = 5)
    pts = Tuple{Float64,Float64}[]
    function flat!(p0, p1, n)
        for k in (isempty(pts) ? 0 : 1):n
            push!(pts, (p0[1] + (p1[1] - p0[1]) * k / n, p0[2] + (p1[2] - p0[2]) * k / n))
        end
    end
    function arc!(c, φ0, φ1, n)
        for k in 1:n
            φ = φ0 + (φ1 - φ0) * k / n
            push!(pts, (c[1] + R * cos(φ), c[2] + R * sin(φ)))
        end
    end
    flat!((B, D - lip), (B, D - R), n_flat)
    arc!((B - R, D - R), 0.0, π / 2, n_corner)
    flat!((B - R, D), (R, D), n_flat)
    arc!((R, D - R), π / 2, π, n_corner)
    flat!((0.0, D - R), (0.0, R), n_flat)
    arc!((R, R), π, 3π / 2, n_corner)
    flat!((R, 0.0), (B - R, 0.0), n_flat)
    arc!((B - R, R), 3π / 2, 2π, n_corner)
    flat!((B, R), (B, lip), n_flat)
    return first.(pts), last.(pts)
end

function torsion_J(section, t; L = 700.0, nz = 56, Cs = TriShellFiniteElement.DEFAULT_SHEAR_RELAXATION,
                   drilling = TriShellFiniteElement.DEFAULT_DRILLING)
    X, Y = section
    G = E / (2 * (1 + ν)); βo = 0.01
    nn = length(X); Z = collect(range(0.0, L, nz + 1))
    Jthin = sum(hypot(X[i+1] - X[i], Y[i+1] - Y[i]) for i in 1:nn-1) * t^3 / 3
    nodes = [Ferrite.Node(Ferrite.Vec(X[i], Y[i], Z[j])) for j in 1:nz+1 for i in 1:nn]
    id(i, j) = (j - 1) * nn + i
    cells = Ferrite.Triangle[]
    for j in 1:nz, i in 1:nn-1
        n1, n2, n3, n4 = id(i, j), id(i, j + 1), id(i + 1, j + 1), id(i + 1, j)
        push!(cells, Ferrite.Triangle((n1, n2, n3))); push!(cells, Ferrite.Triangle((n1, n3, n4)))
    end
    grid = Ferrite.Grid(cells, nodes)
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
    K = TriShellFiniteElement.assemble_global_Ke!(K, dh, qr1, qr3, ip3, ip6, E, ν, t; Cs, drilling)
    K0 = copy(K)
    Xc = sum(X) / nn; Yc = sum(Y) / nn
    end0 = Set(id(i, 1) for i in 1:nn); endL = Set(id(i, nz + 1) for i in 1:nn)
    ch = ConstraintHandler(dh)
    add!(ch, Dirichlet(:u, end0, (x, t_) -> [0.0, 0.0], [1, 2]))
    add!(ch, Dirichlet(:u, Set([id(1, 1)]), (x, t_) -> [0.0], [3]))
    add!(ch, Dirichlet(:u, endL, (x, t_) -> [-βo * (x[2] - Yc), βo * (x[1] - Xc)], [1, 2]))
    close!(ch); update!(ch, 0.0)
    f = zeros(ndofs(dh)); apply!(K, f, ch); u = K \ f; apply!(u, ch)
    R = K0 * u
    T = sum((grid.nodes[n].x[1] - Xc) * R[node_to_dofs[n][2]] - (grid.nodes[n].x[2] - Yc) * R[node_to_dofs[n][1]] for n in endL)
    Fnet = hypot(sum(R[node_to_dofs[n][1]] for n in endL), sum(R[node_to_dofs[n][2]] for n in endL)) / abs(T)
    return (; J = T * L / (G * βo), Jthin, Fnet)
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

        # drilling dof: Hughes–Brezzi (default) vs. the legacy diagonal penalty
        @test TriShellFiniteElement.DEFAULT_DRILLING == :hughes_brezzi
        kp = TriShellFiniteElement.local_elastic_stiffness_matrix!(qr1, qr3, ip3, ip6, E, ν, 1.0, x; Cs = 0.0, drilling = :penalty)
        nd = setdiff(1:18, [6, 12, 18])
        w = [3, 4, 5, 9, 10, 11, 15, 16, 17]
        @test kp[w, w] ≈ ke0[w, w] rtol = 1e-12                            # plate block untouched
        @test count(abs.(eigvals(Symmetric(kp[nd, nd]))) .< 1e-9 * maximum(kp)) == 6   # rigid body modes
        ev = eigvals(Symmetric(ke0))
        @test count(abs.(ev) .< 1e-9 * maximum(ev)) == 6                   # HB: six zero modes incl. drilling
        d = zeros(18)                                                      # rigid in-plane rotation about the origin
        for i in 1:3
            d[6(i-1)+1] = -x[i][2]; d[6(i-1)+2] = x[i][1]; d[6(i-1)+6] = 1.0
        end
        @test norm(ke0 * d) < 1e-9 * maximum(abs, ke0)
        @test_throws ArgumentError TriShellFiniteElement.local_elastic_stiffness_matrix!(qr1, qr3, ip3, ip6, E, ν, 1.0, x; drilling = :none)
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

    @testset "torsion of folded strips (drilling dof treatment)" begin
        # Saint-Venant J from a warping-free static twist (Moen 2008, Sec. 4.2.7.3.2.3): T = G J β′.
        # Thin-walled J = Σ b t³/3; the exact rectangle value carries the free-edge factor 1 − 0.63 t/b.
        tt = 1.0
        flat = torsion_J(strip_section(0.0, 76.2), tt)
        @test flat.J / flat.Jthin ≈ 1 - 0.63 * tt / 76.2 rtol = 0.03
        @test flat.Fnet < 1e-8
        @test torsion_J(strip_section(0.0, 76.2), tt; drilling = :penalty).J ≈ flat.J rtol = 1e-3
        # folding the strip in two at any angle must not change J (the twisting moment crosses the fold)
        for α in (5.0, 30.0, 90.0)
            fold = torsion_J(strip_section(α, 76.2), tt)
            @test fold.J ≈ flat.J rtol = 0.01
            @test fold.Fnet < 1e-8
        end
        # the legacy absolute penalty leaves a spurious net end force (no longer pure torsion)
        @test torsion_J(strip_section(90.0, 76.2), tt; drilling = :penalty).Fnet > 1e-4
        # rounded lipped C: J ≈ Σbt³/3 (exact Saint-Venant is 0.9956 Σbt³/3) with Hughes–Brezzi;
        # the penalty answer is wrong in a mesh- and rotation-center-dependent way (too stiff for
        # quadrilaterals, too soft for this triangle mesh) and carries a net end force
        c = torsion_J(lipped_c_section(), 1.88)
        @test c.J / c.Jthin ≈ 0.99 rtol = 0.03
        @test c.Fnet < 1e-6
        pc = torsion_J(lipped_c_section(), 1.88; drilling = :penalty)
        @test abs(pc.J / c.J - 1) > 0.03
        @test pc.Fnet > 1e-4
    end
end
