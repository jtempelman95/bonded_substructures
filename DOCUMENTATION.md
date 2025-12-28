# Bonded Substructures - Comprehensive Documentation

## Table of Contents
1. [Overview](#overview)
2. [Theoretical Foundation](#theoretical-foundation)
3. [Implementation](#implementation)
4. [Physical Consistency](#physical-consistency)
5. [Contact Enforcement](#contact-enforcement)
6. [Usage Guide](#usage-guide)
7. [Examples](#examples)
8. [Testing](#testing)

---

## Overview

### Purpose
This repository provides a computational framework for analyzing bonded structures with potential interface disbonds using:
- **Mesh generation** with gmsh
- **Craig-Bampton reduced order modeling** (ROM)
- **Contact mechanics** at disbond interfaces
- **Time-domain** and **frequency-domain** response analysis

### Key Features
- ✅ 3D bonded plate geometries with through-thickness material stacking
- ✅ Circular and rectangular disbond regions
- ✅ Craig-Bampton dynamic substructuring for model reduction
- ✅ Contact constraints for disbond interface (penalty method)
- ✅ Harmonic and transient dynamic analysis
- ✅ Visualization with PyVista and Matplotlib

### Coordinate System
**CRITICAL:** The bond interface is in the **x-y plane** (horizontal), with materials stacked in the **z-direction** (through-thickness):

```
z │  ┌─────────────────┐  ← Top surface (z = t1 + t2)
  │  │  Material 2      │
  │  │  (Coating)       │
t1│  ├─────────────────┤  ← Bond interface (x-y plane at z = t1)
  │  │  Material 1      │
  │  │  (Substrate)     │
0 │  └─────────────────┘  ← Bottom surface (z = 0)
      ←──── x, y ────→
```

- **x-direction**: Width (in-plane)
- **y-direction**: Length (in-plane)
- **z-direction**: Through-thickness (stacking)

### Default Dimensions
- **In-plane**: 1 ft × 1 ft (0.3048 m × 0.3048 m)
- **Substrate (t1)**: 0.15 in (3.81 mm)
- **Coating (t2)**: 0.10 in (2.54 mm)
- **Total thickness**: 0.25 in (6.35 mm)
- **Aspect ratio**: 48:1 (valid for thin plate theory)

---

## Theoretical Foundation

### 1. Craig-Bampton Reduction

#### Background
The Craig-Bampton method [1] is a **component mode synthesis** technique that reduces the size of large finite element models while preserving interface physics.

#### Mathematical Formulation

For a structure partitioned into substructures, the displacement of each substructure is:

```
u = [u_i]  ← Interior DOFs
    [u_b]  ← Boundary/Interface DOFs
```

The Craig-Bampton transformation expresses interior DOFs using:
1. **Constraint modes** (Ψ_c): Unit displacement at each boundary DOF with interior fixed
2. **Fixed-interface normal modes** (Φ_n): Free vibration modes with boundary fixed

**Transformation:**
```
u_i = Φ_n * q + Ψ_c * u_b
```

where q are **modal coordinates**.

**Reduced System:**
```
[M_r]{ü} + [K_r]{u} = {F_r}

where:
M_r = T^T * M * T    (Reduced mass matrix)
K_r = T^T * K * T    (Reduced stiffness matrix)

T = [Φ_n  Ψ_c]       (Transformation matrix)
    [ 0    I  ]
```

#### Benefits
- **Massive DOF reduction**: 10,000+ DOFs → 100s DOFs
- **Preserved interface physics**: All boundary DOFs retained
- **Accurate dynamics**: Captures dominant modes
- **Efficient nonlinear analysis**: Small system enables iterative contact solvers

### 2. Contact Mechanics at Disbond Interface

#### Physical Requirements
At a disbond, the interface must satisfy:

1. **Unilateral constraint** (Signorini condition):
   ```
   gap ≥ 0           (No penetration)
   p_n ≤ 0           (Compression only, no tension)
   gap * p_n = 0     (Complementarity)
   ```

2. **Friction** (optional, Coulomb law):
   ```
   |p_t| ≤ μ|p_n|    (Stick condition)
   ```

where:
- gap: Separation between surfaces
- p_n: Normal contact pressure
- p_t: Tangential (friction) stress
- μ: Friction coefficient

#### Current Implementation: Penalty Method

The current implementation uses a **penalty method** (linear approximation):

```
F_contact = k_contact * gap    (if gap < 0)
         = 0                    (if gap ≥ 0)
```

where k_contact = 1×10^10 N/m (very stiff spring).

**Limitations:**
- ⚠️ Linear approximation (actual contact is nonlinear)
- ⚠️ Always applies force (doesn't model true unilateral behavior)
- ⚠️ Requires nonlinear solver for true contact

**Future Extension:**
For true contact enforcement, implement:
```python
while not converged:
    # 1. Compute trial displacement
    u_trial = solve(K_r * u = F)

    # 2. Check contact constraints
    gap = compute_gap(u_trial)

    # 3. Apply contact forces
    if gap < 0:
        F_contact = enforce_no_penetration()
        K_contact = compute_contact_stiffness()
    else:
        F_contact = 0

    # 4. Update and iterate
    u_new = solve((K_r + K_contact) * u = F + F_contact)
```

### 3. Governing Equations

#### Structural Dynamics
```
M * ü + C * u̇ + K * u = F(t)
```

where:
- M: Mass matrix (kg)
- C: Damping matrix (N⋅s/m)
- K: Stiffness matrix (N/m)
- F: Applied force (N)
- u: Displacement (m)

#### Material Constitutive Relations

**Isotropic (e.g., Aluminum):**
```
σ = E/(1-ν²) * [1    ν    0  ] * [ε_x ]
                [ν    1    0  ]   [ε_y ]
                [0    0  (1-ν)/2]  [γ_xy]
```

**Orthotropic (e.g., Carbon/Epoxy):**
```
σ = [E1/(1-ν12*ν21)    ν12*E2/(1-ν12*ν21)    0   ] * [ε_1]
    [ν21*E1/(1-ν12*ν21)    E2/(1-ν12*ν21)    0   ]   [ε_2]
    [0                      0               G12 ]   [γ_12]
```

where:
- E: Young's modulus (Pa)
- ν: Poisson's ratio
- G: Shear modulus (Pa)
- Subscripts 1,2,3: Material principal directions

### 4. Time Integration: Newmark-β Method

For transient analysis:

```
u_{n+1} = u_n + Δt * v_n + (Δt)²/2 * [(1-2β)*a_n + 2β*a_{n+1}]
v_{n+1} = v_n + Δt * [(1-γ)*a_n + γ*a_{n+1}]
```

**Default parameters:**
- β = 0.25 (constant acceleration)
- γ = 0.50 (unconditionally stable)

---

## Implementation

### Repository Structure

```
bonded_substructures/
├── src/bonded_substructures/
│   ├── geometry/
│   │   ├── base.py                   # Abstract base class
│   │   └── rectangle.py              # Bonded plate geometry
│   ├── materials/
│   │   └── properties.py             # Material definitions
│   ├── mesh/
│   │   └── markers.py                # Physical tag management
│   ├── rom/
│   │   └── craig_bampton.py          # Craig-Bampton reduction
│   └── visualization/
│       ├── mesh_plot.py              # Mesh visualization
│       └── results_plot.py           # Results plotting
├── examples/
│   ├── 01_bonded_plate_basic.py      # Basic plate
│   ├── 02_bonded_plate_disbond.py    # Plate with disbond
│   ├── 03_visualize_mesh.py          # Visualization
│   ├── 07_wide_plate_disbond.py      # Large-scale example
│   └── 08_time_domain_response.py    # Transient dynamics
├── tests/
│   └── test_contact_enforcement.py   # Physical consistency tests
└── DOCUMENTATION.md                  # This file
```

### Key Classes

#### `BondedRectangle`
```python
from bonded_substructures.geometry import BondedRectangle

plate = BondedRectangle(
    width=0.3048,           # 1 ft (x-direction)
    length=0.3048,          # 1 ft (y-direction)
    t1=0.00381,             # 0.15 in (substrate thickness)
    t2=0.00254,             # 0.10 in (coating thickness)
    material_1=ALUMINUM_7075_T6,
    material_2=CARBON_EPOXY_UD,
    mesh_size=0.025,        # 25 mm elements
)

# Add disbond
plate.add_disbond(
    position=(0.1524, 0.1524, 0.00381),  # Center at interface
    size=0.0254,                          # 1 in radius
    shape="circular"
)

# Generate mesh
plate.generate_mesh()
plate.save_mesh("bonded_plate.msh")
```

#### `CraigBamptonReduction`
```python
from bonded_substructures.rom.craig_bampton import CraigBamptonReduction

cb = CraigBamptonReduction(
    mesh_file="bonded_plate_disbond.msh",
    material_tags={"Material_1": 1, "Material_2": 2}
)

# Perform reduction
cb.load_mesh()
cb.partition_substructures()
cb.identify_interface_dofs()
cb.assemble_matrices()
cb.compute_modes(n_modes=50)

# Get reduced system
K_r, M_r = cb.assemble_reduced_system()

# Apply contact
cb.apply_contact_constraints(contact_stiffness=1e10)
```

### Physical Tags (gmsh)

| Tag | Description | Dimension |
|-----|-------------|-----------|
| 1 | Material 1 (Substrate) | Volume (3) |
| 2 | Material 2 (Coating) | Volume (3) |
| 30 | Disbond Region | Volume (3) |
| 23 | Top Surface | Surface (2) |

---

## Physical Consistency

### Validated Properties

✅ **Geometry**
- Dimensions: 1 ft × 1 ft × 0.25 in
- Aspect ratio: 48:1 (thin plate valid)
- Interface at z = 0.15 in

✅ **Materials**
- Aluminum 7075-T6: E = 71.7 GPa, ν = 0.33, ρ = 2810 kg/m³
- Carbon/Epoxy UD: E1 = 171 GPa, E2 = 9.1 GPa, ρ = 1570 kg/m³

✅ **Coordinate System**
- x-y: In-plane (0 to 0.3048 m)
- z: Through-thickness (0 to 0.00635 m)
- Bond interface: x-y plane at z = 3.81 mm

### Physical Checks

Run automated tests:
```bash
python tests/test_contact_enforcement.py
```

Tests verify:
1. Geometry consistency (dimensions, aspect ratios)
2. Material property bounds (E > 0, 0 < ν < 0.5)
3. Coordinate system orientation
4. Mesh quality metrics

---

## Contact Enforcement

### Current Implementation

The penalty method applies contact stiffness at interface DOFs:

```python
K_contact = k_penalty * I    # At interface DOFs
```

where k_penalty = 1×10^10 N/m (very stiff).

### Verification Tests

#### Test 1: Stiffness Application
✅ Verifies K_reduced changes after contact application
✅ Checks ||K_with_contact - K_reduced|| > 0

#### Test 2: Matrix Conditioning
✅ Computes condition number: κ(K) = λ_max / λ_min
✅ Ensures κ < 10^15 (well-conditioned)

#### Test 3: Symmetry
✅ Verifies ||K - K^T|| / ||K|| < 10^-10

### Limitations

⚠️ **Current limitations:**
1. **Linear approximation**: True contact is nonlinear (unilateral)
2. **No gap detection**: Penalty always active (not true contact)
3. **No friction**: Only normal contact modeled

⚠️ **Future work needed:**
1. Implement nonlinear contact solver (complementarity problem)
2. Add gap detection and active set algorithm
3. Include Coulomb friction model
4. Implement augmented Lagrangian method

---

## Usage Guide

### Quick Start

```bash
# 1. Install dependencies
conda env create -f environment.yml
conda activate bonded_substructures

# 2. Run basic example
python examples/01_bonded_plate_basic.py

# 3. Run disbond example
python examples/02_bonded_plate_disbond.py

# 4. Visualize meshes
python examples/03_visualize_mesh.py

# 5. Run all core examples
./run_core_examples.sh
```

### Typical Workflow

```python
# Step 1: Create geometry
from bonded_substructures.geometry import BondedRectangle

plate = BondedRectangle(width=0.3048, length=0.3048,
                        t1=0.00381, t2=0.00254)
plate.add_disbond(position=(0.1524, 0.1524, 0.00381),
                  size=0.0254, shape="circular")
plate.generate_mesh()
plate.save_mesh("plate.msh")

# Step 2: Craig-Bampton reduction
from bonded_substructures.rom.craig_bampton import CraigBamptonReduction

cb = CraigBamptonReduction("plate.msh", {"Material_1": 1, "Material_2": 2})
cb.load_mesh()
cb.partition_substructures()
cb.identify_interface_dofs()
cb.assemble_matrices()
cb.compute_modes(n_modes=50)
K_r, M_r = cb.assemble_reduced_system()
cb.apply_contact_constraints(contact_stiffness=1e10)

# Step 3: Time-domain simulation
from examples.08_time_domain_response import newmark_beta

# Define loading
def force(t):
    F = np.zeros(K_r.shape[0])
    if t < 0.001:  # 1 ms impulse
        F[0] = 10000 * (1 - t/0.001)
    return F

# Integrate
t, u, v, a = newmark_beta(M_r, C_r, K_r, force, dt=0.0001, t_end=0.05)

# Step 4: Post-process
import matplotlib.pyplot as plt

plt.plot(t*1000, u[0,:]*1000)
plt.xlabel('Time (ms)')
plt.ylabel('Displacement (mm)')
plt.show()
```

---

## Examples

### Example 01: Basic Plate
- **Purpose**: Generate simple bonded plate mesh
- **Geometry**: 1 ft × 1 ft × 0.25 in, no disbond
- **Output**: bonded_plate_basic.msh (398 nodes, 1929 elements)

### Example 02: Plate with Disbond
- **Purpose**: Add circular disbond at interface
- **Disbond**: 1 in radius (2.2% of interface area)
- **Output**: bonded_plate_disbond.msh (629 nodes, 3811 elements)

### Example 03: Visualization
- **Purpose**: Generate mesh visualizations
- **Methods**: PyVista (3D interactive), Matplotlib (2D plots)
- **Output**: PNG files in examples/output/

### Example 08: Time-Domain Response
- **Purpose**: Transient dynamics with Craig-Bampton ROM
- **Analysis**: Newmark-β time integration
- **Loading**: Impulse (impact), step, harmonic
- **Output**: Displacement, velocity, acceleration histories

---

## Testing

### Physical Consistency Tests

```bash
python tests/test_contact_enforcement.py
```

**Tests performed:**
1. ✅ Geometry consistency (dimensions, total thickness)
2. ✅ Material property bounds
3. ✅ Coordinate system verification
4. ✅ Aspect ratio checks (thin plate theory)
5. ✅ Contact stiffness application
6. ✅ Matrix symmetry and conditioning

**Expected output:**
```
======================================================================
Physical Consistency Tests
======================================================================

Test 1: Geometry Consistency
  ✅ PASS: Geometry is consistent

Test 2: Material Properties
  ✅ PASS: Material properties are physical

Test 3: Coordinate System
  ✅ PASS: Coordinate system is consistent

Test 4: Aspect Ratios
  ✅ Thin plate assumption valid (aspect ratio > 10)
```

### Notes on Craig-Bampton Test

For very thin meshes (t2 << mesh_size), Material 2 may have insufficient interior DOFs for mode extraction. This is **expected behavior** and indicates:
- Mesh refinement needed for thinner coatings
- ROM still valid (constraint modes capture interface physics)
- Contact model still applicable

---

## References

[1] R.R. Craig, M.C.C. Bampton, "Coupling of substructures for dynamic analyses," AIAA Journal, Vol. 6, No. 7, 1968.

[2] K.L. Johnson, "Contact Mechanics," Cambridge University Press, 1985.

[3] P. Wriggers, "Computational Contact Mechanics," 2nd ed., Springer, 2006.

[4] T.J.R. Hughes, "The Finite Element Method: Linear Static and Dynamic Finite Element Analysis," Dover, 2000.

---

## Contact

For questions or issues, please open an issue on GitHub.

**Generated**: 2025-12-28
**Version**: 1.0
**Coordinate System**: x-y plane bond interface, z-direction stacking
