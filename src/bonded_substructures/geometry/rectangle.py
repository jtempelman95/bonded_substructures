"""Rectangle geometry with bonded materials and optional disbond."""

from typing import Optional, Tuple

import gmsh
import numpy as np

from bonded_substructures.geometry.base import BondedGeometry
from bonded_substructures.materials.properties import MaterialProperties
from bonded_substructures.mesh.markers import PhysicalTag


class BondedRectangle(BondedGeometry):
    """2D rectangular domain with two bonded materials.

    The rectangle is divided horizontally into two regions:
    - Bottom region: material_1 (substrate)
    - Top region: material_2 (coating)

    An optional disbond region can be added at the interface.

    Attributes:
        width: Total width of the rectangle
        height_1: Height of material 1 (substrate)
        height_2: Height of material 2 (coating)
        material_1: Properties of substrate material
        material_2: Properties of coating material
        mesh_size: Characteristic mesh element size
    """

    def __init__(
        self,
        width: float,
        height_1: float,
        height_2: float,
        material_1: MaterialProperties,
        material_2: MaterialProperties,
        mesh_size: float = 0.1,
    ):
        """Initialize the bonded rectangle geometry.

        Args:
            width: Total width of the rectangle
            height_1: Height of material 1 (substrate)
            height_2: Height of material 2 (coating)
            material_1: Substrate material properties
            material_2: Coating material properties
            mesh_size: Characteristic mesh element size
        """
        super().__init__(material_1, material_2, mesh_size, dim=2)
        self.width = width
        self.height_1 = height_1
        self.height_2 = height_2
        self.total_height = height_1 + height_2

        # Disbond parameters
        self._has_disbond = False
        self._disbond_position: Optional[Tuple[float, float]] = None
        self._disbond_size: Optional[float] = None
        self._disbond_shape: Optional[str] = None

    def add_disbond(
        self,
        position: Tuple[float, float],
        size: float,
        shape: str = "circular",
    ) -> None:
        """Add a disbond region at the interface.

        The disbond is centered at the specified position on the interface
        between the two materials.

        Args:
            position: (x, y) center position of the disbond
            size: Radius for circular disbond, half-width for rectangular
            shape: Shape of disbond ('circular' or 'rectangular')
        """
        if shape not in ["circular", "rectangular"]:
            raise ValueError("Disbond shape must be 'circular' or 'rectangular'")

        # Validate position is at interface
        if not np.isclose(position[1], self.height_1):
            print(
                f"Warning: Disbond y-position {position[1]} is not at interface "
                f"(y={self.height_1}). Adjusting to interface."
            )
            position = (position[0], self.height_1)

        self._has_disbond = True
        self._disbond_position = position
        self._disbond_size = size
        self._disbond_shape = shape

    def create_geometry(self) -> None:
        """Create the rectangle geometry with bonded materials.

        Creates a 2D rectangle divided into two horizontal regions with
        an optional disbond at the interface.
        """
        # Create points for the rectangle
        # Bottom left, bottom right, top right, top left
        p1 = gmsh.model.geo.addPoint(0, 0, 0, self.mesh_size)
        p2 = gmsh.model.geo.addPoint(self.width, 0, 0, self.mesh_size)
        p3 = gmsh.model.geo.addPoint(self.width, self.total_height, 0, self.mesh_size)
        p4 = gmsh.model.geo.addPoint(0, self.total_height, 0, self.mesh_size)

        # Interface points
        p5 = gmsh.model.geo.addPoint(0, self.height_1, 0, self.mesh_size)
        p6 = gmsh.model.geo.addPoint(self.width, self.height_1, 0, self.mesh_size)

        # Create lines for material 1 (bottom)
        l1 = gmsh.model.geo.addLine(p1, p2)  # Bottom
        l2 = gmsh.model.geo.addLine(p2, p6)  # Right side (partial)
        l3 = gmsh.model.geo.addLine(p6, p5)  # Interface
        l4 = gmsh.model.geo.addLine(p5, p1)  # Left side (partial)

        # Create lines for material 2 (top)
        l5 = gmsh.model.geo.addLine(p5, p6)  # Interface (opposite direction from l3)
        l6 = gmsh.model.geo.addLine(p6, p3)  # Right side (partial)
        l7 = gmsh.model.geo.addLine(p3, p4)  # Top
        l8 = gmsh.model.geo.addLine(p4, p5)  # Left side (partial)

        # Create curve loops and surfaces
        if self._has_disbond:
            # With disbond, we need to create the disbond region and modify the interface
            surface_1, surface_2, disbond_surface = self._create_geometry_with_disbond(
                p1, p2, p3, p4, p5, p6, l1, l2, l3, l4, l5, l6, l7, l8
            )
        else:
            # Without disbond, simple two-surface geometry
            # Bottom surface: p1 -> p2 -> p6 -> p5 -> p1
            curve_loop_1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
            surface_1 = gmsh.model.geo.addPlaneSurface([curve_loop_1])

            # Top surface: p5 -> p6 -> p3 -> p4 -> p5
            curve_loop_2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
            surface_2 = gmsh.model.geo.addPlaneSurface([curve_loop_2])

        gmsh.model.geo.synchronize()

        # Define physical groups
        gmsh.model.addPhysicalGroup(2, [surface_1], PhysicalTag.MATERIAL_1)
        gmsh.model.setPhysicalName(2, PhysicalTag.MATERIAL_1, self.material_1.name)

        gmsh.model.addPhysicalGroup(2, [surface_2], PhysicalTag.MATERIAL_2)
        gmsh.model.setPhysicalName(2, PhysicalTag.MATERIAL_2, self.material_2.name)

        # Boundaries
        gmsh.model.addPhysicalGroup(1, [l1], PhysicalTag.BOUNDARY_BOTTOM)
        gmsh.model.setPhysicalName(1, PhysicalTag.BOUNDARY_BOTTOM, "Bottom")

        gmsh.model.addPhysicalGroup(1, [l7], PhysicalTag.BOUNDARY_TOP)
        gmsh.model.setPhysicalName(1, PhysicalTag.BOUNDARY_TOP, "Top")

        gmsh.model.addPhysicalGroup(1, [l4, l8], PhysicalTag.BOUNDARY_LEFT)
        gmsh.model.setPhysicalName(1, PhysicalTag.BOUNDARY_LEFT, "Left")

        gmsh.model.addPhysicalGroup(1, [l2, l6], PhysicalTag.BOUNDARY_RIGHT)
        gmsh.model.setPhysicalName(1, PhysicalTag.BOUNDARY_RIGHT, "Right")

        # Interface
        if self._has_disbond:
            gmsh.model.addPhysicalGroup(2, [disbond_surface], PhysicalTag.DISBOND_REGION)
            gmsh.model.setPhysicalName(2, PhysicalTag.DISBOND_REGION, "Disbond")
        else:
            gmsh.model.addPhysicalGroup(1, [l3], PhysicalTag.INTERFACE_BONDED)
            gmsh.model.setPhysicalName(1, PhysicalTag.INTERFACE_BONDED, "Bonded Interface")

    def _create_geometry_with_disbond(
        self, p1, p2, p3, p4, p5, p6, l1, l2, l3, l4, l5, l6, l7, l8
    ) -> Tuple[int, int, int]:
        """Create geometry with a disbond region.

        For the initial version, we create a simple embedded disk at the interface
        that will be meshed separately and tagged as the disbond region.

        Args:
            p1-p6: Point tags
            l1-l8: Line tags

        Returns:
            Tuple of (surface_1, surface_2, disbond_surface) tags
        """
        # Create the two main rectangles first (same as without disbond)
        curve_loop_1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
        surface_1 = gmsh.model.geo.addPlaneSurface([curve_loop_1])

        curve_loop_2 = gmsh.model.geo.addCurveLoop([l5, l6, l7, l8])
        surface_2 = gmsh.model.geo.addPlaneSurface([curve_loop_2])

        # Create a thin disbond region embedded at the interface using OCC
        gmsh.model.geo.synchronize()

        x_center, y_center = self._disbond_position
        size = self._disbond_size

        # Create a thin rectangular disbond region at the interface
        # This represents a delamination zone
        disbond_thickness = self.mesh_size * 0.1  # Very thin layer
        y_bottom = y_center - disbond_thickness / 2
        y_top = y_center + disbond_thickness / 2

        if self._disbond_shape == "circular":
            # Create circular disbond using OCC
            disbond_surface = gmsh.model.occ.addDisk(x_center, y_center, 0, size, size)
        else:  # rectangular
            width_disbond = size * 2
            x_corner = x_center - size
            disbond_surface = gmsh.model.occ.addRectangle(
                x_corner, y_bottom, 0, width_disbond, disbond_thickness
            )

        gmsh.model.occ.synchronize()

        # Note: The disbond is created as a separate surface that overlaps the interface
        # In a more sophisticated version, we would use boolean operations to
        # properly split the domains. For now, this creates a marked region.

        return surface_1, surface_2, disbond_surface

    def get_mesh_info(self) -> dict:
        """Get information about the generated mesh.

        Returns:
            Dictionary with mesh statistics
        """
        if not gmsh.isInitialized():
            return {"error": "Mesh not generated"}

        num_nodes = len(gmsh.model.mesh.getNodes()[0])
        num_elements = len(gmsh.model.mesh.getElements()[2][0])

        return {
            "num_nodes": num_nodes,
            "num_elements": num_elements,
            "width": self.width,
            "height_1": self.height_1,
            "height_2": self.height_2,
            "total_height": self.total_height,
            "has_disbond": self._has_disbond,
        }
