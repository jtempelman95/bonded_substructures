# Bonded Substructures - AI Coding Guidelines

## Project Overview
This is a Python package for mesh generation and finite element analysis of bonded materials with disbond regions. It uses gmsh for meshing, FEniCSx for FE analysis, and supports reduced order modeling (ROM) techniques like Craig-Bampton and harmonic balance.

## Architecture
- **geometry/**: Abstract base `BondedGeometry` with concrete implementations (`BondedRectangle`, `BondedCylinder`)
- **materials/**: `MaterialProperties` dataclass supporting isotropic (E, nu) and orthotropic materials
- **mesh/**: Physical tag system and mesh utilities
- **rom/**: Craig-Bampton reduction and harmonic balance analysis
- **solver/**: Linear FE solver
- **visualization/**: Plotting with pyvista/matplotlib

## Key Patterns
### Geometry Creation
Use context managers for automatic gmsh cleanup:
```python
with BondedRectangle(width=10.0, height_1=2.0, height_2=1.0,
                     material_1=ALUMINUM_6061_T6, material_2=CARBON_EPOXY_UD) as geom:
    geom.generate_mesh()
    geom.save_mesh("mesh.msh")
```

### Material Definition
```python
# Isotropic
aluminum = MaterialProperties(name="Aluminum", E=70e9, nu=0.33, rho=2700)

# Orthotropic
composite = MaterialProperties(name="Composite", E1=150e9, E2=10e9, E3=10e9,
                              G12=5e9, G13=5e9, G23=3e9,
                              nu12=0.3, nu13=0.3, nu23=0.4, rho=1600)
```

### Physical Tags
- Materials: 1-9 (MATERIAL_1=1, MATERIAL_2=2)
- Interfaces: 10-19 (INTERFACE_BONDED=10, INTERFACE_DISBOND=11)
- Boundaries: 20-29 (BOUNDARY_LEFT=20, etc.)
- Special: 30-39 (DISBOND_REGION=30)

### Disbond Addition
```python
geom.add_disbond(position=(5.0, 2.0, 1.0), size=1.5, shape="circular")
```

## Development Workflow
- **Install**: `poetry install` (includes dev dependencies)
- **Activate env**: `poetry shell`
- **Run examples**: `python examples/01_rectangle_basic.py`
- **Test**: `pytest` (runs with coverage)
- **Format**: `black .`
- **Lint**: `ruff check .`
- **Build**: `poetry build`

## gmsh Usage
- Initialize with `gmsh.initialize()`
- Use `gmsh.model.occ` for geometry operations
- `gmsh.model.occ.fragment()` for disbond regions
- Always call `gmsh.finalize()` or use context manager
- Physical groups via `gmsh.model.addPhysicalGroup()`

## Testing
- Fixtures in `conftest.py` provide common geometries
- Test mesh files saved to temp directories
- Use pytest-cov for coverage reporting

## Dependencies
- Core: numpy, scipy, gmsh, pyvista, matplotlib
- Optional: fenics-dolfinx (for FE analysis)
- Dev: pytest, black, ruff, ipython

## File Organization
- Examples in `examples/` with descriptive names (01_*, 02_*)
- Tests mirror source structure in `tests/`
- Meshes saved as `.msh` files in workspace root or examples/output/

## Common Pitfalls
- Always finalize gmsh to avoid resource leaks
- Disbond positions must be at material interfaces
- Use absolute paths for file operations
- Check material property validity (isotropic vs orthotropic)