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
        depth: float = 1.0,  # Thickness of the plate in z-direction
    ):
        """Initialize the bonded rectangle geometry.

        Args:
            width: Total width of the rectangle (x-direction)
            height_1: Height of material 1 substrate (y-direction)
            height_2: Height of material 2 coating (y-direction)
            material_1: Substrate material properties
            material_2: Coating material properties
            mesh_size: Characteristic mesh element size
            depth: Thickness/depth of the plate (z-direction)
        """
        super().__init__(material_1, material_2, mesh_size, dim=3)  # Changed to 3D
        self.width = width
        self.height_1 = height_1
        self.height_2 = height_2
        self.total_height = height_1 + height_2
        self.depth = depth

        # Disbond parameters
        self._has_disbond = False
        self._disbond_position: Optional[Tuple[float, float, float]] = None  # Now 3D
        self._disbond_size: Optional[float] = None
        self._disbond_shape: Optional[str] = None

    def add_disbond(
        self,
        position: Tuple[float, float, float],
        size: float,
        shape: str = "circular",
    ) -> None:
        """Add a disbond region at the interface.

        The disbond is centered at the specified position on the interface
        between the two materials.

        Args:
            position: (x, y, z) center position of the disbond
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
            position = (position[0], self.height_1, position[2])

        self._has_disbond = True
        self._disbond_position = position
        self._disbond_size = size
        self._disbond_shape = shape

    def create_geometry(self) -> None:
        """Create 3D plate geometry with bonded materials.

        Creates two 3D boxes (plates) stacked on top of each other with
        an optional disbond volume at the interface.
        """
        # Use OCC kernel for 3D box creation
        gmsh.model.occ.synchronize()

        if self._has_disbond:
            # Create geometry with disbond using boolean operations
            volume_1, volume_2, disbond_volume = self._create_3d_geometry_with_disbond()
        else:
            # Create simple two-box geometry
            # Material 1 (substrate): bottom box
            volume_1 = gmsh.model.occ.addBox(0, 0, 0, self.width, self.height_1, self.depth)

            # Material 2 (coating): top box
            volume_2 = gmsh.model.occ.addBox(
                0, self.height_1, 0, self.width, self.height_2, self.depth
            )

            disbond_volume = None

        gmsh.model.occ.synchronize()

        # Define physical groups for volumes
        gmsh.model.addPhysicalGroup(3, [volume_1], PhysicalTag.MATERIAL_1)
        gmsh.model.setPhysicalName(3, PhysicalTag.MATERIAL_1, self.material_1.name)

        gmsh.model.addPhysicalGroup(3, [volume_2], PhysicalTag.MATERIAL_2)
        gmsh.model.setPhysicalName(3, PhysicalTag.MATERIAL_2, self.material_2.name)

        if disbond_volume is not None:
            gmsh.model.addPhysicalGroup(3, [disbond_volume], PhysicalTag.DISBOND_REGION)
            gmsh.model.setPhysicalName(3, PhysicalTag.DISBOND_REGION, "Disbond")

        # Get boundary surfaces for physical groups
        # This is more complex in 3D - we'll mark key surfaces
        all_surfaces = gmsh.model.occ.getEntities(dim=2)

        # Bottom surface (y=0)
        bottom_surfaces = []
        # Top surface (y=total_height)
        top_surfaces = []

        for dim, tag in all_surfaces:
            com = gmsh.model.occ.getCenterOfMass(dim, tag)
            if np.isclose(com[1], 0, atol=1e-6):  # Bottom
                bottom_surfaces.append(tag)
            elif np.isclose(com[1], self.total_height, atol=1e-6):  # Top
                top_surfaces.append(tag)

        if bottom_surfaces:
            gmsh.model.addPhysicalGroup(2, bottom_surfaces, PhysicalTag.BOUNDARY_BOTTOM)
            gmsh.model.setPhysicalName(2, PhysicalTag.BOUNDARY_BOTTOM, "Bottom")

        if top_surfaces:
            gmsh.model.addPhysicalGroup(2, top_surfaces, PhysicalTag.BOUNDARY_TOP)
            gmsh.model.setPhysicalName(2, PhysicalTag.BOUNDARY_TOP, "Top")

    def _create_3d_geometry_with_disbond(self) -> Tuple[int, int, int]:
        """Create 3D geometry with a disbond volume at the interface.

        Returns:
            Tuple of (volume_1, volume_2, disbond_volume) tags
        """
        x_center, y_center, z_center = self._disbond_position
        size = self._disbond_size

        # Create a thin disbond layer at the interface
        disbond_thickness = self.mesh_size * 0.2  # Thin layer

        # Create the disbond volume
        if self._disbond_shape == "circular":
            # Create cylindrical disbond (disk extruded in z-direction)
            disbond_volume = gmsh.model.occ.addCylinder(
                x_center, y_center - disbond_thickness / 2, z_center,
                0, disbond_thickness, 0,  # Cylinder along y-axis
                size  # radius
            )
        else:  # rectangular
            # Create rectangular box disbond
            width_disbond = size * 2
            depth_disbond = size * 2
            disbond_volume = gmsh.model.occ.addBox(
                x_center - size,
                y_center - disbond_thickness / 2,
                z_center - size,
                width_disbond,
                disbond_thickness,
                depth_disbond
            )

        # Create the two material boxes
        box1 = gmsh.model.occ.addBox(0, 0, 0, self.width, self.height_1, self.depth)
        box2 = gmsh.model.occ.addBox(0, self.height_1, 0, self.width, self.height_2, self.depth)

        # Use fragment to split the boxes at the disbond
        gmsh.model.occ.synchronize()
        out, out_map = gmsh.model.occ.fragment(
            [(3, box1), (3, box2)],
            [(3, disbond_volume)]
        )

        gmsh.model.occ.synchronize()

        # Identify the resulting volumes
        volumes = [tag for dim, tag in out if dim == 3]

        # Classify volumes by their y-centroid
        volume_1 = None
        volume_2 = None
        disbond = None

        for vol in volumes:
            com = gmsh.model.occ.getCenterOfMass(3, vol)
            if com[1] < self.height_1 - disbond_thickness:  # Below interface
                volume_1 = vol
            elif com[1] > self.height_1 + disbond_thickness:  # Above interface
                volume_2 = vol
            else:  # At interface
                disbond = vol

        # Fallback if identification fails
        if volume_1 is None or volume_2 is None or disbond is None:
            print("Warning: Could not properly identify disbond volumes.")
            if len(volumes) >= 3:
                return volumes[0], volumes[1], volumes[2]
            else:
                return box1, box2, disbond_volume

        return volume_1, volume_2, disbond

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
