# TriShellFiniteElement.jl

`TriShellFiniteElement.jl` is a Julia package that implements a **3-node triangular shell finite element** with drilling degrees of freedom (DOF). The element combines membrane, plate bending, and transverse shear stiffness contributions into an 18-DOF shell element suitable for linear static and buckling analysis of shell structures.

---

## Key Features

- 3-node triangular shell element with 6 DOF per node (u, v, w, θₓ, θᵧ, θᵤ)
- Composite element formulation combining:
  - In-plane membrane stiffness
  - Out-of-plane plate bending stiffness (Mindlin-Reissner)
  - Transverse shear stiffness with 5/6 shear correction factor
  - Drilling DOF: Hughes–Brezzi term γ∫(θᵤ − ½(∂v/∂x − ∂u/∂y))² dA (default), or a legacy penalty stiffness of 1/100 of the smallest rotational diagonal term (`drilling = :penalty`)
- Geometric stiffness matrix for linear buckling analysis
- Integration with [Ferrite.jl](https://ferrite-fem.github.io/Ferrite.jl/) for mesh handling and DOF management
- Local-to-global coordinate transformation for arbitrarily oriented shell elements in 3D space

---

## Element Formulation

The element has **18 DOF total** (3 nodes × 6 DOF/node):

```
Node DOF: [u, v, w, θₓ, θᵧ, θᵤ]
```

- **u, v** — in-plane translations (membrane)
- **w** — out-of-plane translation (bending/shear)
- **θₓ, θᵧ** — bending rotations
- **θᵤ** — drilling rotation (in-plane rotation)

The local element stiffness is assembled from:

| Contribution | Matrix size | Notes |
|---|---|---|
| Membrane | 6×6 | Plane stress, 3 nodes × 2 DOF |
| Bending | 9×9 | Plate bending, 3 nodes × 3 DOF (w, θₓ, θᵧ); function returns 18×18 pre-padded with zeros |
| Shear | 18×18 | Transverse shear with drilling DOF, static condensation applied |
| Drilling | 18×18 | Hughes–Brezzi term on the corner functions (default) or 3 diagonal penalty terms |

### Drilling degree of freedom

The drilling rotation θᵤ has no stiffness of its own in a Mindlin shell. `drilling = :hughes_brezzi` (default, `DEFAULT_DRILLING`) adds the Hughes & Brezzi (1989) term γ∫(θᵤ − ½(∂v/∂x − ∂u/∂y))² dA, γ = `drilling_gamma` × G t (`DEFAULT_DRILLING_GAMMA = 1.0`), integrated with the 3-point rule on the linear corner functions; it ties θᵤ to the in-plane rotation of the membrane field, so a rigid in-plane rotation is energy-free and the 18×18 matrix has exactly six zero-energy modes. `drilling = :penalty` is the legacy diagonal penalty of 1/100 of the smallest rotational diagonal term.

The choice does not matter for flat plates but decides the torsional behavior of folded and curved meshes: in Saint-Venant torsion every plate strip rotates in-plane at the rate β′h, which the absolute penalty resists, so the torsion constant of a lipped C comes out far too large and grows with corner refinement; without any drilling stiffness a strip's twisting moment cannot cross a fold and every fold acts as a free edge. The Hughes–Brezzi term transfers the moment across folds and gives J within about 1 % of the exact Saint-Venant value. `test/runtests.jl` ("torsion of folded strips") twists a flat strip, the same strip folded at 5°, 30° and 90°, and a rounded lipped C with a warping-free static twist (Moen 2008, Sec. 4.2.7.3.2.3) and checks these statements. See QuadShellFiniteElement.jl for the same option on the quadrilateral element.

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/runtosolve/TriShellFiniteElement.jl")
```

---

## Basic Usage

### Element-Level Stiffness Matrices

```julia
using Ferrite, Tensors, TriShellFiniteElement

E = 200_000.0   # Young's modulus (MPa)
ν = 0.30        # Poisson's ratio
t = 1.0         # Shell thickness (mm)

# Element node coordinates in local 2D system
x = [Tensors.Vec((0.0, 0.0)),
     Tensors.Vec((100.0, 0.0)),
     Tensors.Vec((0.0, 100.0))]

ip3 = TriShellFiniteElement.IP3()
qr1 = QuadratureRule{RefTriangle}(1)
cv  = CellValues(qr1, ip3, ip3)
reinit!(cv, x)

# Constitutive matrices
Dm = TriShellFiniteElement.calculate_membrane_constitutive_matrix(E, ν, t)
Db = TriShellFiniteElement.calculate_bending_constitutive_matrix(E, ν, t)
Ds = TriShellFiniteElement.calculate_shear_constitutive_matrix(E, ν, t)

# Element stiffness matrices
Ke_m = TriShellFiniteElement.calculate_element_membrane_stiffness_matrix(Dm, cv)  # 6×6
Ke_b = TriShellFiniteElement.calculate_element_bending_stiffness_matrix(Db, cv)   # 18×18
Ke_s = TriShellFiniteElement.calculate_element_shear_stiffness_matrix(Ds, cv)     # 18×18
```

### Global Stiffness Assembly

```julia
using Ferrite, TriShellFiniteElement

E = 200_000.0
ν = 0.30
t = 1.0

# Build mesh and DOF handler
grid = generate_grid(Triangle, (10, 10), Vec((0.0, 0.0, 0.0)), Vec((1000.0, 1000.0, 0.0)))
ip   = Lagrange{RefTriangle, 1}()
dh   = DofHandler(grid)
add!(dh, :u, ip^3)   # u, v, w translations
add!(dh, :θ, ip^3)   # θₓ, θᵧ, θᵤ rotations
close!(dh)

# Set up quadrature and custom interpolations
ip3 = TriShellFiniteElement.IP3()
ip6 = TriShellFiniteElement.IP6()
qr1 = QuadratureRule{RefTriangle}(1)
qr3 = QuadratureRule{RefTriangle}(2)

# Assemble global elastic stiffness
K = allocate_matrix(dh)
K = TriShellFiniteElement.assemble_global_Ke!(K, dh, qr1, qr3, ip3, ip6, E, ν, t)
```

### Buckling Analysis

After solving a linear static problem for stresses σXX, σYY, τXY:

```julia
# Assemble geometric stiffness for buckling
Kg = allocate_matrix(dh)
Kg = TriShellFiniteElement.assemble_global_Kg!(Kg, dh, qr1, ip3, σXX, σYY, τXY)  # σ·t per element, element local frame

# Solve generalized eigenvalue problem: K·φ = λ·Kg·φ
# using your preferred eigensolver
```

---

## Shear Relaxation Factor `Cs`

The element uses Mindlin–Reissner bending with a statically condensed transverse shear stiffness.
For thin shells meshed with elements much larger than the thickness, the unrelaxed shear stiffness
dominates and the element **shear-locks**: bending is too stiff and buckling loads are too high.
To control this the shear stiffness is scaled by `1 / (1 + Cs * alpha)`, where `alpha` is the ratio
of the element shear to bending rotational stiffness (a Tessler–Hughes type relaxation).

`Cs` is a keyword argument of `assemble_global_Ke!` and `local_elastic_stiffness_matrix!` with default
`TriShellFiniteElement.DEFAULT_SHEAR_RELAXATION = 0.2`, calibrated on plate buckling benchmarks
(`test/runtests.jl`). Simply supported square plate in uniform compression, exact k = 4.0,
thickness/width = 1/105:

| Elements across width | `Cs = 0.0` | `Cs = 0.2` |
|---|---|---|
| 8  | 4.97 | 4.09 |
| 16 | 4.08 | 4.00 |
| 24 | 4.01 | 3.98 |

With `Cs = 0` the result at a fixed mesh also depends strongly on thickness (k = 10.6 at t/b = 1/460
on an 8×8 mesh, 5.0 at 16×16); with `Cs = 0.2` it is thickness independent (k = 4.10, 4.02, 4.00 at
8, 16, 24 across for t/b = 1/460 and 4.03, 3.94, 3.92 for t/b = 1/46).

> **Bug fixed in v0.5.0.** Up to v0.4.0 the ratio `alpha` was computed over dofs 4:9 of the
> condensed 9-dof plate block, which mixed the `w` translations into the rotational stiffness
> sums. Since the two have different units, `alpha` and hence the relaxation depended on the unit
> system: the calibration above held in millimetres, but the same plate in inches gave k = 2.93
> at 16 elements across. `alpha` is now the ratio on the rotational dofs only (as in
> QuadShellFiniteElement.jl) and the results are unit independent; in millimetre units they are
> within 0.3 % of the old ones. Pass `Cs = 0.0` to recover the
unrelaxed element:

```julia
K = TriShellFiniteElement.assemble_global_Ke!(K, dh, qr1, qr3, ip3, ip6, E, ν, t; Cs = 0.0)
```

---

## Tests

```julia
using Pkg; Pkg.test("TriShellFiniteElement")
```

`test/runtests.jl` builds simply supported and clamped rectangular plates in uniform compression,
solves the pre-buckling membrane problem, assembles the geometric stiffness and checks the buckling
coefficient against the classical plate solutions (k = 4.0 simply supported, k ≈ 10.07 clamped),
including thickness independence and invariance to the plate's orientation in 3D. The other files in
`test/` are development scripts and are not run by the test suite.

---

## Elastic–plastic, geometrically nonlinear analysis

Version 0.5.0 adds von Mises plasticity with isotropic hardening and a geometrically nonlinear
internal force vector / tangent stiffness, with the same API in TriShellFiniteElement.jl and
QuadShellFiniteElement.jl (`src/plasticity.jl` and `src/solver.jl` are shared verbatim).

- **Material**: `J2PlaneStress(E, ν, hardening)` with `IsotropicHardening(εp, σ)`, a piecewise-linear
  true stress – equivalent plastic strain table (first point = initial yield, constant beyond the last
  point, like ABAQUS `*PLASTIC`); `IsotropicHardening(σy)` is perfectly plastic and
  `IsotropicHardening(σy, H)` linear hardening. The plane-stress return mapping is the projected
  algorithm of Simo & Taylor (1986) with a safeguarded scalar Newton on the plastic multiplier and
  the algorithmically consistent tangent (checked against finite differences in the tests).
- **Section**: `ShellSection(mat, t; nlayers = 5, rule = :gauss | :simpson)` integrates the layer
  stresses through the thickness into membrane forces and moments and the coupled section tangent;
  `rule = :simpson` with an odd `nlayers` includes the surfaces and the midplane (ABAQUS shell default).
- **Element**: total Lagrangian in each element's initial local frame. The membrane strain is the
  Green–Lagrange strain of the midsurface with the gradients of all three displacements
  (`εx = u,x + ½(u,x² + v,x² + w,x²)`, …), the same quadratic terms the geometric stiffness of the
  linear buckling analysis is built from; the curvatures are linear in the rotations and the
  transverse shear and drilling terms stay linear elastic (with the same `Cs` relaxation and drilling
  treatment as the linear element). The tangent at zero displacement is exactly the matrix of
  `assemble_global_Ke!`. The formulation is meant for moderate rotations (plate and
  thin-walled member collapse); it is not a corotational large-rotation shell.
- **Model and solver**: `PlasticShellModel(dh, section; Cs, drilling, drilling_gamma)` stores the
  element data and plastic history, `assemble_internal!(K, f, model, u)` assembles the tangent and
  internal force, `commit!` / `restore!` accept or discard the state of a step, `element_results`
  returns membrane forces, moments, equivalent plastic strain and layer stresses, and
  `displaced_grid` imposes an imperfection. `solve_nonlinear(model, ch; ...)` is an incremental
  Newton–Raphson driver for a Ferrite `ConstraintHandler` whose Dirichlet values depend on the
  pseudo-time (displacement control) and/or an external load `f_ext(t)`, with a tangent predictor
  at each step, a backtracking line search, adaptive step cutting/growing and optional artificial
  damping (`damping = :auto`, the `*STATIC, STABILIZE` idea with a unit-density lumped mass
  from `lumped_mass`) for mode jumps and snaps under displacement control.

```julia
using TriShellFiniteElement, Ferrite
M = TriShellFiniteElement

hard = IsotropicHardening([0.0, 0.0016, 0.0259, 0.07, 0.1407, 0.1988],      # true plastic strain
                          [33.0, 34.31, 42.14, 51.35, 60.86, 65.5])          # true stress, ksi
mat = J2PlaneStress(29500.0, 0.3, hard)
sec = ShellSection(mat, 0.0346; nlayers = 5, rule = :simpson)
model = PlasticShellModel(dh, sec)              # dh: fields :u, :θ on the (imperfect) mesh

ch = ConstraintHandler(dh)                      # displacement control: edges move ±t·δmax
add!(ch, Dirichlet(:u, perimeter, (x, t) -> [0.0], [3]))
add!(ch, Dirichlet(:u, edge_x0, (x, t) -> [t * δmax], [1]))
add!(ch, Dirichlet(:u, edge_xL, (x, t) -> [-t * δmax], [1]))
close!(ch)

hist = solve_nonlinear(model, ch; reaction_dofs = u_dofs_of_edge_xL, Δt = 0.01, Δt_max = 0.02)
P = -hist.P                                     # reaction (compression positive) vs hist.t
res = element_results(model, hist.u)            # N, M, ε̄p, layer stresses per integration point
```

### Validation: stiffened element loaded to collapse (Moen 2008, Section 7.1)

The simply supported plate of Moen (2008), "Direct Strength Design of Cold-Formed Steel Members
with Perforations", PhD dissertation, Johns Hopkins University, Section 7.1 (h = 3.4 in, L = 27.2 in,
t = 0.0346 in, E = 29500 ksi, ν = 0.3, the true stress–strain curve of Fig. 7.1, loaded edges
displaced uniformly, imperfection = fundamental buckling mode scaled to d₁/t) was run with this
element on a 128 × 16 mesh (4096 elements) with 5 Simpson layers and compared with the ABAQUS S9R5
results of the dissertation. See the study in
`RunToSolve/Projects/OneRack/base_plate_modeling_2025/Tri_Quad_isotropic_hardening_plasticity_CDM_dissertation_plates_2026`
for the scripts, all cases (imperfection series, Fig. 7.2 boundary conditions, slotted hole) and figures.

| quantity | Moen (2008) | TriShellFiniteElement |
|---|---|---|
| elastic buckling Pcr/Pᵧ | 0.33 (Fig. 7.16) | 0.336 |
| peak load, d₁/t = 0.34, P/Pᵧ | 0.50 (Fig. 7.17) | 0.512 |
| δ/t at the peak | ≈ 0.50 | 0.52 |
| effective width at the peak, midplane S11, hₑ/h | 0.51 (0.48–0.55, Fig. 7.21) | 0.55 (0.53–0.58) |

---

## Dependencies

| Package | Role |
|---|---|
| [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) | Mesh, DOF handler, cell values, quadrature |
| [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl) | Coordinate vectors and tensor operations |
| LinearAlgebra | Matrix operations (standard library) |
| SparseArrays | Artificial damping diagonal in the nonlinear solver (standard library) |
| [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) | Small fixed-size vectors and matrices in the plasticity routines |

---
N