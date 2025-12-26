"""Base class for bonded geometry definitions."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import gmsh
import numpy as np

from bonded_substructures.materials.properties import MaterialProperties


class BondedGeometry(ABC):
    """Abstract base class for bonded material geometries.

    This class defines the interface for creating meshes of bonded structures
    with optional disbond regions. Subclasses implement specific geometries
    (rectangle, cylinder, etc.).

    Attributes:
        material_1: Properties of the first material (e.g., substrate)
        material_2: Properties of the second material (e.g., coating)
        mesh_size: Characteristic mesh element size
        dim: Geometric dimension (2 for 2D, 3 for 3D)
    """

    def __init__(
        self,
        material_1: MaterialProperties,
        material_2: MaterialProperties,
        mesh_size: float = 0.1,
        dim: int = 2,
    ):
        """Initialize the bonded geometry.

        Args:
            material_1: First material properties
            material_2: Second material properties
            mesh_size: Characteristic mesh element size
            dim: Geometric dimension (2 or 3)
        """
        self.material_1 = material_1
        self.material_2 = material_2
        self.mesh_size = mesh_size
        self.dim = dim
        self._model_name = self.__class__.__name__

    @abstractmethod
    def create_geometry(self) -> None:
        """Create the geometry using gmsh API.

        This method should:
        1. Initialize gmsh if needed
        2. Create geometric entities (points, lines, surfaces, volumes)
        3. Define physical groups for materials and boundaries
        4. NOT generate the mesh (that's done in generate_mesh)
        """
        pass

    @abstractmethod
    def add_disbond(
        self,
        position: Tuple[float, ...],
        size: float,
        shape: str = "circular",
    ) -> None:
        """Add a disbond region to the geometry.

        Args:
            position: Center position of the disbond (x, y) or (x, y, z)
            size: Characteristic size (radius for circular, side length for rectangular)
            shape: Shape of disbond ('circular' or 'rectangular')
        """
        pass

    def generate_mesh(self) -> None:
        """Generate the finite element mesh.

        This method:
        1. Calls create_geometry() if not already done
        2. Sets mesh size parameters
        3. Generates the mesh
        """
        gmsh.initialize()
        gmsh.model.add(self._model_name)

        self.create_geometry()

        # Set mesh size
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", self.mesh_size * 0.5)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", self.mesh_size * 2.0)

        # Generate mesh
        gmsh.model.mesh.generate(self.dim)

    def save_mesh(self, filename: str) -> None:
        """Save the mesh to a file.

        Args:
            filename: Output filename (e.g., 'mesh.msh')
        """
        output_path = Path(filename)
        gmsh.write(str(output_path))
        print(f"Mesh saved to {output_path.absolute()}")

    def get_physical_tags(self) -> Dict[str, int]:
        """Get dictionary mapping physical group names to tags.

        Returns:
            Dictionary with physical group names as keys and tags as values
        """
        physical_groups = {}
        dim_tags = gmsh.model.getPhysicalGroups()

        for dim, tag in dim_tags:
            name = gmsh.model.getPhysicalName(dim, tag)
            physical_groups[name] = tag

        return physical_groups

    def visualize(self) -> None:
        """Launch gmsh GUI to visualize the geometry/mesh."""
        if gmsh.isInitialized():
            gmsh.fltk.run()
        else:
            print("Gmsh not initialized. Call generate_mesh() first.")

    def finalize(self) -> None:
        """Clean up and finalize gmsh."""
        if gmsh.isInitialized():
            gmsh.finalize()

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit - ensures gmsh is finalized."""
        self.finalize()
