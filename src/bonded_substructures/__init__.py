"""
Bonded Substructures - Mesh generation and analysis for bonded materials.

This package provides tools for creating finite element meshes of bonded materials
(e.g., aluminum cylinder with composite overwrap) with optional disbond regions.
It supports mesh generation with gmsh, FE analysis with dolfinx, and visualization
with pyvista and matplotlib.
"""

__version__ = "0.1.0"

from bonded_substructures.materials.properties import MaterialProperties

__all__ = ["MaterialProperties"]
