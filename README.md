# Bonded Substructures

Mesh generation and analysis tools for bonded materials with optional disbond regions using gmsh and FEniCSx.

## Overview

This package provides tools for creating finite element meshes of bonded materials (e.g., aluminum cylinder with composite overwrap) with optional disbond regions. It supports:

- **Mesh generation** with gmsh for complex bonded geometries
- **Material modeling** for isotropic and orthotropic materials
- **Disbond regions** at material interfaces
- **Visualization** with pyvista and matplotlib
- **Future support** for FE analysis with dolfinx and reduced order modeling via Craig-Bampton

## Features

### Current Features (v0.1.0)

- ✅ 2D rectangle geometry with two bonded materials
- ✅ Optional circular or rectangular disbond regions
- ✅ Predefined material properties (aluminum, steel, carbon/epoxy, glass/epoxy)
- ✅ Mesh generation with gmsh
- ✅ Interactive visualization with pyvista
- ✅ 2D plotting with matplotlib
- ✅ Physical tag system for materials, interfaces, and boundaries
- ✅ Full test suite with pytest

### Planned Features

- 🔄 3D cylindrical geometry for COPV (Composite Overwrapped Pressure Vessel) modeling
- 🔄 FE analysis with dolfinx (linear elasticity)
- 🔄 Craig-Bampton reduced order modeling (ROM)
- 🔄 Nonlinear contact mechanics for disbond analysis
- 🔄 Multi-layer composite laminates with stacking sequences

## Installation

### Prerequisites

- Python 3.9 or higher
- Poetry (recommended) or pip

### Install with Poetry (Recommended)

```bash
# Clone the repository
git clone <repository-url>
cd bonded_substructures

# Install dependencies with Poetry
poetry install

# Activate the virtual environment
poetry shell
```

### Install with pip

```bash
# Clone the repository
git clone <repository-url>
cd bonded_substructures

# Install in development mode
pip install -e .
```

### Dependencies

Core dependencies:
- `gmsh` - Mesh generation
- `fenics-dolfinx` - Finite element analysis (future use)
- `pyvista` - Interactive 3D visualization
- `matplotlib` - 2D plotting
- `numpy` - Numerical operations
- `scipy` - Scientific computing

Development dependencies:
- `pytest` - Testing framework
- `black` - Code formatting
- `ruff` - Linting

## Quick Start

### Example 1: Basic Bonded Rectangle

```python
from bonded_substructures.geometry import BondedRectangle
from bonded_substructures.materials.properties import ALUMINUM_6061_T6, CARBON_EPOXY_UD
from bonded_substructures.visualization import print_mesh_info

# Create geometry
with BondedRectangle(
    width=10.0,
    height_1=2.0,  # Aluminum substrate
    height_2=1.0,  # Carbon/epoxy coating
    material_1=ALUMINUM_6061_T6,
    material_2=CARBON_EPOXY_UD,
    mesh_size=0.5,
) as geom:
    # Generate mesh
    geom.generate_mesh()

    # Print mesh information
    print_mesh_info()

    # Save mesh
    geom.save_mesh("my_mesh.msh")
```

### Example 2: Adding a Disbond Region

```python
from bonded_substructures.geometry import BondedRectangle
from bonded_substructures.materials.properties import ALUMINUM_7075_T6, CARBON_EPOXY_UD

# Create geometry
geom = BondedRectangle(
    width=10.0,
    height_1=2.0,
    height_2=1.0,
    material_1=ALUMINUM_7075_T6,
    material_2=CARBON_EPOXY_UD,
    mesh_size=0.3,
)

# Add circular disbond at interface
geom.add_disbond(
    position=(5.0, 2.0),  # Center at interface
    size=1.5,             # Radius = 1.5 m
    shape="circular"
)

# Generate and save mesh
geom.generate_mesh()
geom.save_mesh("disbond_mesh.msh")
geom.finalize()
```

### Example 3: Visualization

```python
from bonded_substructures.visualization import plot_mesh_matplotlib, plot_mesh_pyvista
import matplotlib.pyplot as plt

# Matplotlib visualization (2D)
fig = plot_mesh_matplotlib("my_mesh.msh", show_physical_groups=True)
plt.show()

# PyVista visualization (interactive 3D)
plotter = plot_mesh_pyvista("my_mesh.msh", show_edges=True)
plotter.show()
```

## Running Examples

The `examples/` directory contains ready-to-run examples:

```bash
# Basic rectangle mesh
python examples/01_rectangle_basic.py

# Rectangle with disbond
python examples/02_rectangle_disbond.py

# Visualize meshes
python examples/03_visualize_mesh.py
```

## Project Structure

```
bonded_substructures/
├── src/bonded_substructures/
│   ├── geometry/         # Mesh generation modules
│   │   ├── base.py       # Abstract base class
│   │   └── rectangle.py  # Rectangle geometry
│   ├── materials/        # Material property definitions
│   ├── mesh/            # Mesh utilities and markers
│   ├── visualization/   # Plotting tools
│   ├── solver/          # FE solver (future)
│   └── rom/             # Reduced order modeling (future)
├── examples/            # Example scripts
├── tests/              # Test suite
└── pyproject.toml      # Project configuration
```

## Material Properties

Predefined materials are available in `bonded_substructures.materials.properties`:

**Isotropic Materials:**
- `ALUMINUM_6061_T6` - Aluminum 6061-T6
- `ALUMINUM_7075_T6` - Aluminum 7075-T6
- `STEEL_4130` - Steel 4130

**Orthotropic Materials:**
- `CARBON_EPOXY_UD` - Unidirectional carbon/epoxy (IM7/8552)
- `GLASS_EPOXY_UD` - Unidirectional glass/epoxy

### Custom Materials

```python
from bonded_substructures.materials.properties import MaterialProperties

# Isotropic material
my_material = MaterialProperties(
    name="Custom Aluminum",
    E=70e9,      # Young's modulus (Pa)
    nu=0.33,     # Poisson's ratio
    rho=2700,    # Density (kg/m³)
)

# Orthotropic material
my_composite = MaterialProperties(
    name="Custom Composite",
    E1=150e9, E2=10e9, E3=10e9,           # Young's moduli (Pa)
    G12=5e9, G13=5e9, G23=3e9,            # Shear moduli (Pa)
    nu12=0.3, nu13=0.3, nu23=0.4,         # Poisson's ratios
    rho=1600,                              # Density (kg/m³)
)
```

## Physical Tags

The package uses a standardized physical tag system for mesh regions:

| Tag | Region | Description |
|-----|--------|-------------|
| 1-9 | Materials | Material regions (MATERIAL_1, MATERIAL_2, ...) |
| 10-19 | Interfaces | Bonded and disbond interfaces |
| 20-29 | Boundaries | External boundaries (left, right, top, bottom, etc.) |
| 30-39 | Special | Disbond regions |

Access tags via `bonded_substructures.mesh.markers.PhysicalTag`.

## Testing

Run the test suite with pytest:

```bash
# Run all tests
pytest

# Run with coverage
pytest --cov=bonded_substructures --cov-report=html

# Run specific test file
pytest tests/test_geometry/test_rectangle.py
```

## Future Development: Craig-Bampton ROM

The architecture is designed to support Craig-Bampton reduced order modeling for efficient nonlinear disbond analysis:

1. **Substructuring**: Partition mesh at disbond interface (substrate + coating)
2. **Reduction**: Apply Craig-Bampton to each substructure
   - Compute constraint modes (interface DOFs)
   - Compute fixed-interface normal modes (interior DOFs)
3. **Coupling**: Preserve interface DOFs for nonlinear contact
4. **Analysis**: Apply contact conditions in reduced space

This approach dramatically reduces DOFs while preserving interface physics for disbond modeling.

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass: `pytest`
5. Format code: `black .`
6. Lint code: `ruff check .`
7. Submit a pull request

## License

[Specify your license here]

## Citation

If you use this software in your research, please cite:

```
[Add citation information]
```

## Contact

[Add contact information]

## Acknowledgments

This package uses:
- [gmsh](https://gmsh.info/) for mesh generation
- [FEniCSx](https://fenicsproject.org/) for finite element analysis
- [PyVista](https://pyvista.org/) for visualization
- [NumPy](https://numpy.org/) and [SciPy](https://scipy.org/) for numerical computing
