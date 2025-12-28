"""Rectangle geometry with bonded materials and optional disbond.

COORDINATE SYSTEM:
- x-y plane: In-plane dimensions (width × length)
- z-direction: Through-thickness (material stacking)
- Bond interface: x-y plane at z = t1
- Material 1 (substrate): 0 ≤ z ≤ t1
- Material 2 (coating): t1 ≤ z ≤ t1 + t2
"""

from typing import Optional, Tuple

import gmsh
import numpy as np

from bonded_substructures.geometry.base import BondedGeometry
from bonded_substructures.materials.properties import MaterialProperties
from bonded_substructures.mesh.markers import PhysicalTag


class BondedRectangle(BondedGeometry):
    """Bonded plate geometry with interface in x-y plane.

    Creates a rectangular bonded plate with two materials stacked
    in the z-direction (through-thickness). The bond interface is
    a horizontal plane at z = t1.

    Physical setup:
        - Material 1 (substrate): Bottom layer, 0 ≤ z ≤ t1
        - Material 2 (coating): Top layer, t1 ≤ z ≤ t1 + t2
        - Bond interface: x-y plane at z = t1
        - Disbond: Region in x-y plane at interface

    Default dimensions (easily modifiable):
        - width: 0.3048 m (1 ft, x-direction)
        - length: 0.3048 m (1 ft, y-direction)
        - t1: 0.00381 m (0.15 in substrate thickness)
        - t2: 0.00254 m (0.10 in coating thickness)
        - total: 0.00635 m (0.25 in total thickness)

    Attributes:
        width: Plate width in x-direction (m)
        length: Plate length in y-direction (m)
        t1: Substrate thickness in z-direction (m)
        t2: Coating thickness in z-direction (m)
        material_1: Properties of substrate material
        material_2: Properties of coating material
        mesh_size: Characteristic mesh element size (m)
    """

    def __init__(
        self,
        width: float = 0.3048,  # 1 ft
        length: float = 0.3048,  # 1 ft
        t1: float = 0.00381,  # 0.15 in substrate
        t2: float = 0.00254,  # 0.10 in coating
        material_1: MaterialProperties = None,
        material_2: MaterialProperties = None,
        mesh_size: float = 0.025,  # 25 mm (about 1 inch)
    ):
        """Initialize the bonded plate geometry.

        Args:
            width: Plate width in x-direction (m), default 0.3048 m (1 ft)
            length: Plate length in y-direction (m), default 0.3048 m (1 ft)
            t1: Substrate thickness in z-direction (m), default 0.00381 m (0.15 in)
            t2: Coating thickness in z-direction (m), default 0.00254 m (0.10 in)
            material_1: Substrate material properties
            material_2: Coating material properties
            mesh_size: Characteristic mesh element size (m), default 0.025 m (25 mm)
        """
        # Default materials if not provided
        if material_1 is None:
            from bonded_substructures.materials import ALUMINUM_7075_T6
            material_1 = ALUMINUM_7075_T6
        if material_2 is None:
            from bonded_substructures.materials import CARBON_EPOXY_UD
            material_2 = CARBON_EPOXY_UD

        super().__init__(material_1, material_2, mesh_size, dim=3)

        self.width = width
        self.length = length
        self.t1 = t1
        self.t2 = t2
        self.total_thickness = t1 + t2

        # Disbond parameters
        self._has_disbond = False
        self._disbond_position: Optional[Tuple[float, float, float]] = None
        self._disbond_size: Optional[float] = None
        self._disbond_shape: Optional[str] = None

    def add_disbond(
        self,
        position: Tuple[float, float, float],
        size: float,
        shape: str = "circular",
    ) -> None:
        """Add a disbond region at the bond interface.

        The disbond is a region in the x-y plane at z = t1 (bond interface).

        Args:
            position: (x, y, z) center position of the disbond
                     Note: z will be automatically set to t1 (bond interface)
            size: Radius for circular disbond, half-width for rectangular (m)
            shape: Shape of disbond ('circular' or 'rectangular')
        """
        if shape not in ["circular", "rectangular"]:
            raise ValueError("Disbond shape must be 'circular' or 'rectangular'")

        # Force disbond to be at interface (z = t1)
        if not np.isclose(position[2], self.t1):
            print(
                f"Warning: Disbond z-position {position[2]} is not at interface "
                f"(z={self.t1}). Adjusting to interface."
            )
            position = (position[0], position[1], self.t1)

        self._has_disbond = True
        self._disbond_position = position
        self._disbond_size = size
        self._disbond_shape = shape

    def create_geometry(self) -> None:
        """Create 3D bonded plate geometry.

        Creates two 3D boxes (plates) stacked in z-direction:
        - Box 1: 0 ≤ x ≤ width, 0 ≤ y ≤ length, 0 ≤ z ≤ t1
        - Box 2: 0 ≤ x ≤ width, 0 ≤ y ≤ length, t1 ≤ z ≤ t1+t2

        If disbond present, uses boolean fragment to split the interface.
        """
        gmsh.model.occ.synchronize()

        if self._has_disbond:
            volume_1, volume_2, disbond_volume = self._create_geometry_with_disbond()
        else:
            # Simple stacked boxes
            volume_1 = gmsh.model.occ.addBox(0, 0, 0, self.width, self.length, self.t1)
            volume_2 = gmsh.model.occ.addBox(0, 0, self.t1, self.width, self.length, self.t2)
            disbond_volume = None

        gmsh.model.occ.synchronize()

        # Define physical groups for volumes (dimension 3)
        gmsh.model.addPhysicalGroup(3, [volume_1], PhysicalTag.MATERIAL_1)
        gmsh.model.setPhysicalName(3, PhysicalTag.MATERIAL_1, self.material_1.name)

        gmsh.model.addPhysicalGroup(3, [volume_2], PhysicalTag.MATERIAL_2)
        gmsh.model.setPhysicalName(3, PhysicalTag.MATERIAL_2, self.material_2.name)

        if disbond_volume is not None:
            gmsh.model.addPhysicalGroup(3, [disbond_volume], PhysicalTag.DISBOND_REGION)
            gmsh.model.setPhysicalName(3, PhysicalTag.DISBOND_REGION, "Disbond")

        # Add boundary surfaces
        self._add_boundary_surfaces()

    def _create_geometry_with_disbond(self) -> Tuple[int, int, int]:
        """Create geometry with disbond region at interface.

        Returns:
            Tuple of (volume_1_tag, volume_2_tag, disbond_volume_tag)
        """
        x_center, y_center, z_center = self._disbond_position
        size = self._disbond_size

        # Disbond is a thin disk in x-y plane at z = t1
        disbond_thickness = self.mesh_size * 0.2  # Thin layer

        if self._disbond_shape == "circular":
            # Create a thin cylinder in x-y plane (axis in z-direction)
            disbond_volume = gmsh.model.occ.addCylinder(
                x_center, y_center, z_center - disbond_thickness / 2,
                0, 0, disbond_thickness,  # Cylinder axis in z-direction
                size  # radius
            )
        else:  # rectangular
            # Create a thin box
            disbond_volume = gmsh.model.occ.addBox(
                x_center - size, y_center - size, z_center - disbond_thickness / 2,
                2 * size, 2 * size, disbond_thickness
            )

        # Create main material volumes
        box1 = gmsh.model.occ.addBox(0, 0, 0, self.width, self.length, self.t1)
        box2 = gmsh.model.occ.addBox(0, 0, self.t1, self.width, self.length, self.t2)

        gmsh.model.occ.synchronize()

        # Use fragment to split materials at disbond interface
        # This creates separate volumes for bonded and disbonded regions
        out, out_map = gmsh.model.occ.fragment(
            [(3, box1), (3, box2)],
            [(3, disbond_volume)]
        )

        gmsh.model.occ.synchronize()

        # Identify volumes by their z-coordinates
        volume_1 = None
        volume_2 = None
        disbond_vol = None

        for dim, tag in out:
            if dim == 3:  # Volume
                # Get bounding box to identify volume
                bbox = gmsh.model.occ.getBoundingBox(dim, tag)
                z_min, z_max = bbox[2], bbox[5]
                z_center = (z_min + z_max) / 2

                # Check volume centroid to classify
                mass = gmsh.model.occ.getMass(dim, tag)
                com = gmsh.model.occ.getCenterOfMass(dim, tag)

                if z_center < self.t1 * 0.5:
                    # Bottom material (substrate)
                    if abs(mass - self.width * self.length * self.t1) < 1e-6:
                        volume_1 = tag  # Main substrate
                elif z_center > self.t1 + self.t2 * 0.5:
                    # Top material (coating)
                    if abs(mass - self.width * self.length * self.t2) < 1e-6:
                        volume_2 = tag  # Main coating
                else:
                    # Near interface - likely disbond
                    if abs(z_center - self.t1) < disbond_thickness:
                        disbond_vol = tag

        # If classification failed, use default ordering
        if volume_1 is None or volume_2 is None:
            print("Warning: Volume classification failed, using default ordering")
            volumes = [tag for dim, tag in out if dim == 3]
            volume_1 = volumes[0] if len(volumes) > 0 else box1
            volume_2 = volumes[1] if len(volumes) > 1 else box2
            disbond_vol = volumes[2] if len(volumes) > 2 else disbond_volume

        return volume_1, volume_2, disbond_vol

    def _add_boundary_surfaces(self) -> None:
        """Add physical groups for boundary surfaces."""
        # Get all surfaces
        surfaces = gmsh.model.occ.getEntities(dim=2)

        # Classify surfaces by position
        bottom_surfaces = []  # z = 0
        top_surfaces = []     # z = t1 + t2
        side_surfaces = []    # x or y boundaries

        for dim, tag in surfaces:
            bbox = gmsh.model.occ.getBoundingBox(dim, tag)
            x_min, y_min, z_min = bbox[0], bbox[1], bbox[2]
            x_max, y_max, z_max = bbox[3], bbox[4], bbox[5]

            # Bottom surface (z = 0)
            if np.isclose(z_min, 0.0) and np.isclose(z_max, 0.0):
                bottom_surfaces.append(tag)
            # Top surface (z = t1 + t2)
            elif np.isclose(z_min, self.total_thickness) and np.isclose(z_max, self.total_thickness):
                top_surfaces.append(tag)

        # Create physical groups for boundaries
        if bottom_surfaces:
            gmsh.model.addPhysicalGroup(2, bottom_surfaces, PhysicalTag.BOUNDARY_BOTTOM)
            gmsh.model.setPhysicalName(2, PhysicalTag.BOUNDARY_BOTTOM, "Bottom")

        if top_surfaces:
            gmsh.model.addPhysicalGroup(2, top_surfaces, PhysicalTag.BOUNDARY_TOP)
            gmsh.model.setPhysicalName(2, PhysicalTag.BOUNDARY_TOP, "Top")
