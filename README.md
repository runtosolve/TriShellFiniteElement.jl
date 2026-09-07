# TriShellFiniteElement.jl

`TriShellFiniteElement.jl` is a Julia package that implements a **3-node triangular shell finite element** with drilling degrees of freedom (DOF). The element combines membrane, plate bending, and transverse shear stiffness contributions into an 18-DOF shell element suitable for linear static and buckling analysis of shell structures.

---

## Key Features

- 3-node triangular shell element with 6 DOF per node (u, v, w, θₓ, θᵧ, θᵤ)
- Composite element formulation combining:
  - In-plane membrane stiffness
  - Out-of-plane plate bending stiffness (Mindlin-Reissner)
  - Transverse shear stiffness with 5/6 shear correction factor
  - Drilling DOF stabilization
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
| 8  | 4.97 | 4.10 |
| 10 | 4.44 | 4.05 |
| 16 | 4.08 | 4.00 |
| 24 | 4.01 | 3.98 |
| 32 | 3.99 | 3.97 |

With `Cs = 0` the result at a fixed mesh also depends strongly on thickness (k = 8.3 at t/b = 1/460
on a 10×10 mesh); with `Cs = 0.2` it is thickness independent. Pass `Cs = 0.0` to recover the
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

## Dependencies

| Package | Role |
|---|---|
| [Ferrite.jl](https://github.com/Ferrite-FEM/Ferrite.jl) | Mesh, DOF handler, cell values, quadrature |
| [Tensors.jl](https://github.com/Ferrite-FEM/Tensors.jl) | Coordinate vectors and tensor operations |
| LinearAlgebra | Matrix operations (standard library) |

---
N