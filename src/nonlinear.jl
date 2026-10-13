# --------------------------------------------------------------------------------------------
# Geometrically nonlinear, elastic–plastic triangular shell element.
#
# Total Lagrangian formulation in each element's initial local frame (the frame of the undeformed,
# possibly imperfect, mesh): the membrane strain is the Green–Lagrange strain of the midsurface
#   εx  = u,x + ½ (u,x² + v,x² + w,x²)
#   εy  = v,y + ½ (u,y² + v,y² + w,y²)
#   γxy = u,y + v,x + (u,x u,y + v,x v,y + w,x w,y)
# (gradients of all three displacement components, as in the geometric stiffness matrix of
# Visy & Ádány 2017 / Moen & Ádány 2025), the curvatures are linear in the rotations
# (κx = θy,x, κy = −θx,y, κxy = θy,y − θx,x) and the transverse shear strains and the drilling
# term are linear elastic. The membrane and bending stresses are coupled through the layered
# elastic–plastic section (plasticity.jl). The consistent tangent is the material tangent plus the
# geometric stiffness of the membrane forces N; at zero displacement it is exactly the linear
# element stiffness of `local_elastic_stiffness_matrix!` (same shear relaxation `Cs` and drilling
# treatment).
#
# Element dof order in the local frame: [u v w θx θy θz] per node (18 dofs).
# --------------------------------------------------------------------------------------------

const IND_FIELD_18 = [1, 2, 3, 7, 8, 9, 13, 14, 15, 4, 5, 6, 10, 11, 12, 16, 17, 18]   # block → Ferrite field order

"""
    TriShellElementData

Per-element data of a [`PlasticShellModel`](@ref): initial local frame, constant shape function
gradients, the linear elastic part of the stiffness (relaxed transverse shear + drilling) and the
converged and pending layer states of the single integration point.
"""
mutable struct TriShellElementData
    T::Matrix{Float64}              # 3×3 local frame
    Te::Matrix{Float64}             # 18×18 block rotation
    dNdx::Vector{Vec{2,Float64}}    # constant gradients of the three linear functions
    area::Float64
    k_el::Matrix{Float64}           # 18×18 linear elastic part: transverse shear (condensed, relaxed) + drilling
    states::Vector{LayerState}      # converged layer states
    pending::Vector{LayerState}     # states of the latest evaluation
end

"""
    PlasticShellModel(dh, section::ShellSection; Cs, drilling, drilling_gamma)

Geometrically nonlinear elastic–plastic shell model on the mesh of the DofHandler `dh` (fields
`:u` and `:θ` on linear triangles, 3D node coordinates; an initial imperfection is simply part of
the node coordinates). `section` holds the material, thickness and through-thickness layers.
`Cs`, `drilling` and `drilling_gamma` have the same meaning as in [`assemble_global_Ke!`](@ref).
Use [`assemble_internal!`](@ref) for the internal force vector and tangent stiffness,
[`solve_nonlinear`](@ref) to run an incremental analysis, and [`element_results`](@ref) for
stresses.
"""
struct PlasticShellModel
    dh::DofHandler
    section::ShellSection
    elements::Vector{TriShellElementData}
    Cs::Float64
    drilling::Symbol
    drilling_gamma::Float64
end

function PlasticShellModel(dh::DofHandler, section::ShellSection; Cs = DEFAULT_SHEAR_RELAXATION,
                           drilling = DEFAULT_DRILLING, drilling_gamma = DEFAULT_DRILLING_GAMMA)
    E = section.mat.E; ν = section.mat.ν; t = section.t
    ip3 = IP3(); ip6 = IP6()
    qr1 = QuadratureRule{RefTriangle}(1); qr3 = QuadratureRule{RefTriangle}(2)
    cv = CellValues(qr1, ip3, ip3)
    nl = nlayers(section)
    elements = Vector{TriShellElementData}(undef, getncells(dh.grid))
    for cell in CellIterator(dh)
        x_global = getcoordinates(cell)
        T = calculation_rotation_matrix(x_global)
        x_local = global_nodal_coords_to_planar_coords(x_global, T)
        reinit!(cv, x_local)
        dNdx = [shape_gradient(cv, 1, i) for i in 1:3]
        area = getdetJdV(cv, 1)

        ke_m, ke_b, ke_s = local_elastic_pieces(qr1, qr3, ip3, ip6, E, ν, t, x_local; Cs)
        induv = [1 2 6 7 11 12]; indwt = [3 4 5 8 9 10 13 14 15]
        ind18 = [1:5; 7:11; 13:17]
        k15 = zeros(15, 15)                     # full linear matrix, only used by the :penalty drilling option
        k15[induv, induv] = ke_m
        k15[indwt, indwt] = ke_s + ke_b
        # linear elastic part of the nonlinear element: transverse shear + drilling (membrane and
        # bending come from the layered section)
        k15_s = zeros(15, 15)
        k15_s[indwt, indwt] = ke_s
        k_el = zeros(18, 18)
        k_el[ind18, ind18] = k15_s
        add_drilling_stiffness!(k_el, k15, qr3, ip3, x_local, E, ν, t, drilling, drilling_gamma)

        Te = rotation_matrix_for_element_stiffness_drilling(T)
        elements[cellid(cell)] = TriShellElementData(T, Te, dNdx, area, k_el,
                                                     [LayerState() for _ in 1:nl], [LayerState() for _ in 1:nl])
    end
    return PlasticShellModel(dh, section, elements, Cs, drilling, drilling_gamma)
end

"""
    tri_internal_force!(fe, ke, el, sec, d)

Internal force `fe` (18) and tangent `ke` (18×18) of one element in its local frame for the local
nodal displacements `d` (block order `[u v w θx θy θz]` per node). Updates `el.pending`.
"""
function tri_internal_force!(fe, ke, el::TriShellElementData, sec::ShellSection, d::AbstractVector)
    dN = el.dNdx
    ux = uy = vx = vy = wx = wy = 0.0
    κx = κy = κxy = 0.0
    @inbounds for i in 1:3
        o = 6 * (i - 1)
        gx = dN[i][1]; gy = dN[i][2]
        ux += gx * d[o+1]; uy += gy * d[o+1]
        vx += gx * d[o+2]; vy += gy * d[o+2]
        wx += gx * d[o+3]; wy += gy * d[o+3]
        κx += gx * d[o+5]
        κy -= gy * d[o+4]
        κxy += gy * d[o+5] - gx * d[o+4]
    end
    εm = SVector(ux + 0.5 * (ux^2 + vx^2 + wx^2),
                 vy + 0.5 * (uy^2 + vy^2 + wy^2),
                 uy + vx + ux * uy + vx * vy + wx * wy)
    κ = SVector(κx, κy, κxy)

    N, M, Dmm, Dmb, Dbb = section_response(sec, εm, κ, el.states, el.pending)

    Bm = _tri_Bm(dN, ux, uy, vx, vy, wx, wy)
    Bb = _tri_Bb(dN)

    A = el.area
    k = A * (Bm' * Dmm * Bm + Bm' * Dmb * Bb + Bb' * Dmb' * Bm + Bb' * Dbb * Bb)
    fl = A * (Bm' * N + Bb' * M)
    @inbounds for j in 1:18, i in 1:18
        ke[i, j] = k[i, j] + el.k_el[i, j]
    end
    @inbounds for i in 1:18
        fe[i] = fl[i]
    end
    mul!(fe, el.k_el, d, 1.0, 1.0)

    # geometric stiffness of the membrane forces on the three translations
    @inbounds for j in 1:3, i in 1:3
        gxi = dN[i][1]; gyi = dN[i][2]; gxj = dN[j][1]; gyj = dN[j][2]
        g = A * (N[1] * gxi * gxj + N[2] * gyi * gyj + N[3] * (gxi * gyj + gyi * gxj))
        for c in 1:3
            ke[6*(i-1)+c, 6*(j-1)+c] += g
        end
    end
    return fe, ke
end

function _tri_Bm(dN, ux, uy, vx, vy, wx, wy)
    B = zeros(MMatrix{3,18,Float64})
    @inbounds for i in 1:3
        o = 6 * (i - 1)
        gx = dN[i][1]; gy = dN[i][2]
        B[1, o+1] = gx * (1 + ux); B[1, o+2] = gx * vx;       B[1, o+3] = gx * wx
        B[2, o+1] = gy * uy;       B[2, o+2] = gy * (1 + vy); B[2, o+3] = gy * wy
        B[3, o+1] = gy * (1 + ux) + gx * uy
        B[3, o+2] = gx * (1 + vy) + gy * vx
        B[3, o+3] = gy * wx + gx * wy
    end
    return SMatrix(B)
end

function _tri_Bb(dN)
    B = zeros(MMatrix{3,18,Float64})
    @inbounds for i in 1:3
        o = 6 * (i - 1)
        gx = dN[i][1]; gy = dN[i][2]
        B[1, o+5] = gx
        B[2, o+4] = -gy
        B[3, o+4] = -gx; B[3, o+5] = gy
    end
    return SMatrix(B)
end

"""
    assemble_internal!(K, f, model::PlasticShellModel, u)

Assemble the global internal force vector `f` and tangent stiffness `K` for the global solution
vector `u` (Ferrite dof order of `model.dh`). The layer states of the latest evaluation are kept
as pending until [`commit!`](@ref) is called.
"""
function assemble_internal!(K, f, model::PlasticShellModel, u::AbstractVector)
    assembler = start_assemble(K, f)
    fe = zeros(18); ke = zeros(18, 18)
    ug = zeros(18); ul = zeros(18)
    fg = zeros(18); kg = zeros(18, 18); tmp = zeros(18, 18)
    fe_f = zeros(18); ke_f = zeros(18, 18)
    for cell in CellIterator(model.dh)
        el = model.elements[cellid(cell)]
        cd = celldofs(cell)
        @inbounds for k in 1:18
            ug[IND_FIELD_18[k]] = u[cd[k]]      # Ferrite field order → nodal block order
        end
        mul!(ul, el.Te', ug)
        tri_internal_force!(fe, ke, el, model.section, ul)
        mul!(fg, el.Te, fe)
        mul!(tmp, ke, el.Te')
        mul!(kg, el.Te, tmp)
        @inbounds for j in 1:18
            fe_f[j] = fg[IND_FIELD_18[j]]
            for i in 1:18
                ke_f[i, j] = kg[IND_FIELD_18[i], IND_FIELD_18[j]]
            end
        end
        assemble!(assembler, cd, ke_f, fe_f)
    end
    return K, f
end

"""
    commit!(model)

Accept the pending layer states of the latest [`assemble_internal!`](@ref) as converged.
"""
function commit!(model::PlasticShellModel)
    for el in model.elements
        copyto!(el.states, el.pending)
    end
    return model
end

"""
    restore!(model)

Discard the pending states (nothing else to restore for the triangle, whose internal dofs are
linear elastic and condensed exactly).
"""
function restore!(model::PlasticShellModel)
    for el in model.elements
        copyto!(el.pending, el.states)
    end
    return model
end

"""
    element_results(model::PlasticShellModel, u)

Section results at the integration point of every element for the solution `u`: a vector of
named tuples `(x, N, M, εm, κ, ε̄p, σ_layers)` with the global coordinates of the point, the
membrane forces and moments (force/length, moment/length) in the element local frame, the section
strains, the maximum equivalent plastic strain over the layers and the layer stresses (ordered as
`model.section.z`). The pending states are not modified.
"""
function element_results(model::PlasticShellModel, u::AbstractVector)
    sec = model.section
    out = Vector{NamedTuple}(undef, getncells(model.dh.grid))
    ug = zeros(18); d = zeros(18)
    scratch = [LayerState() for _ in 1:nlayers(sec)]
    for cell in CellIterator(model.dh)
        el = model.elements[cellid(cell)]
        cd = celldofs(cell)
        for k in 1:18
            ug[IND_FIELD_18[k]] = u[cd[k]]
        end
        mul!(d, el.Te', ug)
        dN = el.dNdx
        ux = uy = vx = vy = wx = wy = 0.0; κx = κy = κxy = 0.0
        for i in 1:3
            o = 6 * (i - 1); gx = dN[i][1]; gy = dN[i][2]
            ux += gx * d[o+1]; uy += gy * d[o+1]; vx += gx * d[o+2]; vy += gy * d[o+2]
            wx += gx * d[o+3]; wy += gy * d[o+3]
            κx += gx * d[o+5]; κy -= gy * d[o+4]; κxy += gy * d[o+5] - gx * d[o+4]
        end
        εm = SVector(ux + 0.5 * (ux^2 + vx^2 + wx^2), vy + 0.5 * (uy^2 + vy^2 + wy^2),
                     uy + vx + ux * uy + vx * vy + wx * wy)
        κ = SVector(κx, κy, κxy)
        N, M, _, _, _ = section_response(sec, εm, κ, el.states, scratch)
        σl = layer_stresses(sec, εm, κ, el.states)
        xg = getcoordinates(cell)
        x = sum(xg) / 3
        out[cellid(cell)] = (; x, N, M, εm, κ, ε̄p = maximum(s.ε̄p for s in el.states), σ_layers = σl)
    end
    return out
end

"""
    lumped_mass(model::PlasticShellModel)

Lumped mass vector (Ferrite dof order) of the shell with unit density: each element's mass
`t × area` is split equally between its nodes on the three translations, with the rotary inertia
`t² / 12` of a plate per unit mass on the three rotations. Used by the artificial damping of
[`solve_nonlinear`](@ref).
"""
function lumped_mass(model::PlasticShellModel)
    m = zeros(ndofs(model.dh))
    t = model.section.t
    for cell in CellIterator(model.dh)
        el = model.elements[cellid(cell)]
        me = t * el.area / 3
        cd = celldofs(cell)
        for i in 1:3, k in 1:3
            m[cd[3(i-1)+k]] += me
            m[cd[9+3(i-1)+k]] += me * t^2 / 12
        end
    end
    return m
end

"""
    displaced_grid(grid, Δx)

Copy of `grid` with every node `i` moved by `Δx[i]` (a `Vec{3}`), e.g. to impose an initial
geometric imperfection from a scaled buckling mode.
"""
function displaced_grid(grid::Grid, Δx::AbstractVector)
    nodes = [Ferrite.Node(n.x + Δx[i]) for (i, n) in enumerate(grid.nodes)]
    return Grid(copy(grid.cells), nodes; cellsets = copy(Ferrite.getcellsets(grid)),
                nodesets = copy(Ferrite.getnodesets(grid)), facetsets = copy(Ferrite.getfacetsets(grid)))
end
