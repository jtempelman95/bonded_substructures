"""Hollow cylinder geometry with radial bonding interface.

This module implements a COPV (Composite Overwrapped Pressure Vessel) style
geometry with concentric cylindrical shells bonded at a radial interface.
"""

from typing import Optional, Tuple
import numpy as np
import gmsh

from .base import BondedGeometry
from ..materials import MaterialProperties, ALUMINUM_7075_T6, CARBON_EPOXY_UD
from ..mesh.markers import PhysicalTag


class BondedCylinder(BondedGeometry):
    """Hollow cylinder geometry with radial bonding interface (COPV-style).

    Physical setup:
        - Material 1 (substrate): Inner cylindrical shell
          Radial extent: r_inner ≤ r ≤ r_inner + t1
        - Material 2 (coating): Outer cylindrical shell
          Radial extent: r_inner + t1 ≤ r ≤ r_inner + t1 + t2
        - Bond interface: Cylindrical surface at r = r_inner + t1
        - Disbond: Optional patch region on cylindrical interface
        - Coordinate system: Cylinder axis along z-direction, origin at bottom center

    Default dimensions (COPV scale):
        - radius_inner: 0.45 m (~1.5 ft) - inner bore radius
        - t1: 0.003 m (3 mm, ~0.12 in) - aluminum liner thickness
        - t2: 0.020 m (20 mm, ~0.79 in) - composite overwrap thickness
        - height: 1.8 m (~6 ft) - axial length
        - Total outer radius: 0.473 m (~1.55 ft)

    Attributes:
        radius_inner (float): Inner bore radius (m)
        t1 (float): Substrate thickness in radial direction (m)
        t2 (float): Coating thickness in radial direction (m)
        height (float): Axial height (m)
        radius_interface (float): Bond interface radius = radius_inner + t1 (m)
        radius_outer (float): Outer surface radius = radius_inner + t1 + t2 (m)
    """

    def __init__(
        self,
        radius_inner: float = 0.45,
        t1: float = 0.003,
        t2: float = 0.020,
        height: float = 1.8,
        material_1: Optional[MaterialProperties] = None,
        material_2: Optional[MaterialProperties] = None,
        mesh_size: float = 0.05,
    ):
        """Initialize hollow cylinder geometry.

        Args:
            radius_inner: Inner bore radius (m). Default 0.45 m (~1.5 ft)
            t1: Substrate thickness in radial direction (m). Default 0.003 m (3 mm)
            t2: Coating thickness in radial direction (m). Default 0.020 m (20 mm)
            height: Axial height (m). Default 1.8 m (~6 ft)
            material_1: Substrate material. Default ALUMINUM_7075_T6
            material_2: Coating material. Default CARBON_EPOXY_UD
            mesh_size: Characteristic mesh element size (m). Default 0.05 m (50 mm)
        """
        # Set default materials if not provided
        if material_1 is None:
            material_1 = ALUMINUM_7075_T6
        if material_2 is None:
            material_2 = CARBON_EPOXY_UD

        # Initialize base class (3D geometry)
        super().__init__(material_1, material_2, mesh_size, dim=3)

        # Store geometry parameters
        self.radius_inner = radius_inner
        self.t1 = t1
        self.t2 = t2
        self.height = height

        # Compute derived parameters
        self.radius_interface = radius_inner + t1  # Bond interface radius
        self.radius_outer = radius_inner + t1 + t2  # Outer surface radius

        # Disbond parameters (initially no disbond)
        self._has_disbond = False
        self._disbond_position: Optional[Tuple[float, float, float]] = None
        self._disbond_size: Optional[Tuple[float, float]] = None
        self._disbond_shape: Optional[str] = None

    def create_geometry(self) -> None:
        """Create hollow cylinder geometry with gmsh.

        Creates concentric cylindrical shells for substrate and coating materials.
        If disbond is present, calls _create_geometry_with_disbond().
        Otherwise creates simple bonded cylinders.
        """
        gmsh.model.occ.synchronize()

        if self._has_disbond:
            vol_1, vol_2, disbond_vol = self._create_geometry_with_disbond()
        else:
            # Create simple bonded cylinders
            vol_1, vol_2 = self._create_simple_cylinders()
            disbond_vol = None

        gmsh.model.occ.synchronize()

        # Add physical groups for materials (dimension 3 = volumes)
        gmsh.model.addPhysicalGroup(3, [vol_1], PhysicalTag.MATERIAL_1)
        gmsh.model.setPhysicalName(3, PhysicalTag.MATERIAL_1, self.material_1.name)

        gmsh.model.addPhysicalGroup(3, [vol_2], PhysicalTag.MATERIAL_2)
        gmsh.model.setPhysicalName(3, PhysicalTag.MATERIAL_2, self.material_2.name)

        if disbond_vol is not None:
            gmsh.model.addPhysicalGroup(3, [disbond_vol], PhysicalTag.DISBOND_REGION)
            gmsh.model.setPhysicalName(3, PhysicalTag.DISBOND_REGION, "Disbond")

        # Add boundary surfaces
        self._add_boundary_surfaces()

    def _create_simple_cylinders(self) -> Tuple[int, int]:
        """Create simple bonded cylinders without disbond.

        Returns:
            Tuple of (substrate_volume_tag, coating_volume_tag)
        """
        # Create three concentric cylinders:
        # 1. Inner bore (to be subtracted)
        # 2. Substrate outer surface (interface)
        # 3. Coating outer surface

        cylinder_bore = gmsh.model.occ.addCylinder(
            0, 0, 0,               # Center at origin
            0, 0, self.height,     # Axis along z-direction
            self.radius_inner      # Inner bore radius
        )

        cylinder_interface = gmsh.model.occ.addCylinder(
            0, 0, 0,
            0, 0, self.height,
            self.radius_interface  # Substrate outer = interface
        )

        cylinder_outer = gmsh.model.occ.addCylinder(
            0, 0, 0,
            0, 0, self.height,
            self.radius_outer      # Coating outer surface
        )

        gmsh.model.occ.synchronize()

        # Create substrate: interface cylinder - bore
        out_substrate, _ = gmsh.model.occ.cut(
            [(3, cylinder_interface)],
            [(3, cylinder_bore)],
            removeObject=True,
            removeTool=True
        )
        vol_substrate = out_substrate[0][1]

        # Create coating: outer cylinder - interface cylinder
        out_coating, _ = gmsh.model.occ.cut(
            [(3, cylinder_outer)],
            [(3, cylinder_interface)],
            removeObject=True,
            removeTool=False  # Keep interface for bonding
        )
        vol_coating = out_coating[0][1]

        return (vol_substrate, vol_coating)

    def _add_boundary_surfaces(self) -> None:
        """Add physical groups for cylinder boundary surfaces.

        Classifies surfaces by radial and axial position:
        - Inner surface: r ≈ radius_inner (inner bore)
        - Outer surface: r ≈ radius_outer (outer surface)
        - Bottom: z ≈ 0 (bottom end cap)
        - Top: z ≈ height (top end cap)
        """
        # Get all surface entities
        surfaces = gmsh.model.occ.getEntities(dim=2)

        # Lists to collect surface tags by boundary type
        inner_surfaces = []    # Inner bore
        outer_surfaces = []    # Outer surface
        bottom_surfaces = []   # Bottom end cap
        top_surfaces = []      # Top end cap

        for dim, tag in surfaces:
            # Get bounding box
            bbox = gmsh.model.occ.getBoundingBox(dim, tag)
            x_min, y_min, z_min = bbox[0], bbox[1], bbox[2]
            x_max, y_max, z_max = bbox[3], bbox[4], bbox[5]

            # Compute radial extents (distance from z-axis)
            r_min = np.sqrt(x_min**2 + y_min**2)
            r_max = np.sqrt(x_max**2 + y_max**2)
            r_avg = (r_min + r_max) / 2

            # Compute axial position
            z_avg = (z_min + z_max) / 2

            # Classify surface by position
            # Use relative tolerance for radial, absolute for axial
            if np.isclose(r_avg, self.radius_inner, rtol=0.05):
                # Inner bore surface
                inner_surfaces.append(tag)
            elif np.isclose(r_avg, self.radius_outer, rtol=0.05):
                # Outer surface
                outer_surfaces.append(tag)
            elif np.isclose(z_avg, 0.0, atol=self.mesh_size):
                # Bottom end cap
                bottom_surfaces.append(tag)
            elif np.isclose(z_avg, self.height, atol=self.mesh_size):
                # Top end cap
                top_surfaces.append(tag)

        # Create physical groups for boundaries (dimension 2 = surfaces)
        if inner_surfaces:
            gmsh.model.addPhysicalGroup(2, inner_surfaces, PhysicalTag.BOUNDARY_INNER)
            gmsh.model.setPhysicalName(2, PhysicalTag.BOUNDARY_INNER, "Inner")

        if outer_surfaces:
            gmsh.model.addPhysicalGroup(2, outer_surfaces, PhysicalTag.BOUNDARY_OUTER)
            gmsh.model.setPhysicalName(2, PhysicalTag.BOUNDARY_OUTER, "Outer")

        if bottom_surfaces:
            gmsh.model.addPhysicalGroup(2, bottom_surfaces, PhysicalTag.BOUNDARY_END_1)
            gmsh.model.setPhysicalName(2, PhysicalTag.BOUNDARY_END_1, "Bottom")

        if top_surfaces:
            gmsh.model.addPhysicalGroup(2, top_surfaces, PhysicalTag.BOUNDARY_END_2)
            gmsh.model.setPhysicalName(2, PhysicalTag.BOUNDARY_END_2, "Top")

    def add_disbond(
        self,
        position: Tuple[float, float, float],
        size: Tuple[float, float],
        shape: str = "rectangular",
    ) -> None:
        """Add disbond patch on cylindrical interface.

        The disbond is a patch region on the cylindrical bonding interface.
        Position is specified in cylindrical coordinates (theta, z).

        Args:
            position: (theta_center_deg, z_center, _) where:
                - theta_center_deg: Circumferential angle in degrees (0-360)
                - z_center: Axial position (0 to height)
                - Third component ignored (placeholder)
            size: (theta_extent_deg, z_extent) where:
                - theta_extent_deg: Angular extent in degrees
                - z_extent: Axial extent (height of patch)
            shape: Shape of disbond. Only "rectangular" currently supported.

        Raises:
            ValueError: If shape is not "rectangular" or z-position is invalid
        """
        if shape != "rectangular":
            raise ValueError("Cylinder disbond currently only supports 'rectangular' shape")

        theta_center_deg, z_center, _ = position
        theta_extent_deg, z_extent = size

        # Validate axial position is within cylinder height
        if not (0 <= z_center <= self.height):
            raise ValueError(
                f"Disbond z-position {z_center:.3f} m is outside cylinder height "
                f"[0, {self.height:.3f}] m"
            )

        # Validate angular extent
        if not (0 < theta_extent_deg <= 360):
            raise ValueError(
                f"Disbond angular extent {theta_extent_deg} deg must be in range (0, 360]"
            )

        # Validate axial extent
        if not (0 < z_extent <= self.height):
            raise ValueError(
                f"Disbond axial extent {z_extent:.3f} m must be in range (0, {self.height:.3f}]"
            )

        # Store disbond parameters
        self._has_disbond = True
        self._disbond_position = (theta_center_deg, z_center, 0)
        self._disbond_size = (theta_extent_deg, z_extent)
        self._disbond_shape = shape

    def _create_geometry_with_disbond(self) -> Tuple[int, int, int]:
        """Create cylinder geometry with disbond patch at interface.

        Creates a thin annular cylindrical section at the interface representing
        the disbond region. For simplicity, uses full circumferential band.
        Angular sector cutting is noted as future enhancement.

        Returns:
            Tuple of (substrate_volume, coating_volume, disbond_volume)
        """
        theta_center_deg, z_center, _ = self._disbond_position
        theta_extent_deg, z_extent = self._disbond_size

        # Note: Angular sector cutting not implemented yet
        # Using full circumferential band for now
        if theta_extent_deg < 360:
            print(
                f"Note: Angular sector disbonds not yet implemented. "
                f"Creating full circumferential disbond band instead."
            )

        # Create thin cylindrical shell at interface for disbond
        disbond_thickness = self.mesh_size * 0.2  # Thin layer
        r_disbond_inner = self.radius_interface - disbond_thickness / 2
        r_disbond_outer = self.radius_interface + disbond_thickness / 2

        # Create full annular section
        disbond_outer_cyl = gmsh.model.occ.addCylinder(
            0, 0, z_center - z_extent / 2,
            0, 0, z_extent,
            r_disbond_outer
        )

        disbond_inner_cyl = gmsh.model.occ.addCylinder(
            0, 0, z_center - z_extent / 2,
            0, 0, z_extent,
            r_disbond_inner
        )

        gmsh.model.occ.synchronize()

        # Create annular shell
        out_annulus, _ = gmsh.model.occ.cut(
            [(3, disbond_outer_cyl)],
            [(3, disbond_inner_cyl)],
            removeObject=True,
            removeTool=True
        )
        disbond_annulus = out_annulus[0][1]

        gmsh.model.occ.synchronize()

        # Create substrate and coating cylinders separately to avoid fragment issues
        # Substrate: Inner shell
        cyl_substrate_outer = gmsh.model.occ.addCylinder(
            0, 0, 0, 0, 0, self.height, self.radius_interface
        )
        cyl_bore = gmsh.model.occ.addCylinder(
            0, 0, 0, 0, 0, self.height, self.radius_inner
        )

        gmsh.model.occ.synchronize()

        # Cut bore from substrate
        out_substrate, _ = gmsh.model.occ.cut(
            [(3, cyl_substrate_outer)],
            [(3, cyl_bore)],
            removeObject=True,
            removeTool=True
        )
        vol_substrate = out_substrate[0][1]

        # Coating: Outer shell
        cyl_coating_outer = gmsh.model.occ.addCylinder(
            0, 0, 0, 0, 0, self.height, self.radius_outer
        )
        cyl_coating_inner = gmsh.model.occ.addCylinder(
            0, 0, 0, 0, 0, self.height, self.radius_interface
        )

        gmsh.model.occ.synchronize()

        # Cut inner from outer
        out_coating, _ = gmsh.model.occ.cut(
            [(3, cyl_coating_outer)],
            [(3, cyl_coating_inner)],
            removeObject=True,
            removeTool=True
        )
        vol_coating = out_coating[0][1]

        vol_disbond = disbond_annulus

        gmsh.model.occ.synchronize()

        return (vol_substrate, vol_coating, vol_disbond)

    def _cut_to_angular_sector(
        self, annulus_tag: int, theta_center_deg: float, theta_extent_deg: float
    ) -> int:
        """Cut annular shell to angular sector using cutting planes.

        Args:
            annulus_tag: gmsh tag of annular cylinder to cut
            theta_center_deg: Center angle in degrees
            theta_extent_deg: Angular extent in degrees

        Returns:
            Tag of resulting sector volume
        """
        # Convert to radians
        theta_center = np.deg2rad(theta_center_deg)
        theta_half = np.deg2rad(theta_extent_deg / 2)

        # Compute plane normals for sector boundaries
        # Plane 1: at theta_center - theta_half
        theta_1 = theta_center - theta_half
        normal_1 = np.array([-np.sin(theta_1), np.cos(theta_1), 0])

        # Plane 2: at theta_center + theta_half
        theta_2 = theta_center + theta_half
        normal_2 = np.array([np.sin(theta_2), -np.cos(theta_2), 0])

        # Create large boxes oriented with planes for cutting
        # (gmsh OCC doesn't have direct plane cutting, use box approximation)
        box_size = 2 * max(self.radius_outer, self.height)

        # For simplicity, if angular extent is large (> 180deg), keep full annulus
        # Otherwise this becomes complex with plane intersections
        if theta_extent_deg >= 180:
            return annulus_tag

        # Simplified approach: Keep full annulus for now
        # Full angular sector cutting would require complex plane intersection logic
        # This is noted as future enhancement
        print(
            f"Warning: Angular sector cutting not fully implemented. "
            f"Using full {theta_extent_deg:.1f}° circumferential disbond."
        )

        return annulus_tag

    def _classify_volumes_with_disbond(self, volumes) -> Tuple[int, int, int]:
        """Classify volumes after fragment operation with disbond.

        Args:
            volumes: List of (dim, tag) tuples from fragment operation

        Returns:
            Tuple of (substrate_tag, coating_tag, disbond_tag)
        """
        # Extract volume tags
        vol_tags = [tag for dim, tag in volumes if dim == 3]

        # Classify by computing volume and radial position
        vol_substrate = None
        vol_coating = None
        vol_disbond = None

        for tag in vol_tags:
            # Get volume properties
            mass = gmsh.model.occ.getMass(3, tag)
            com = gmsh.model.occ.getCenterOfMass(3, tag)

            # Compute average radius
            r_avg = np.sqrt(com[0]**2 + com[1]**2)

            # Classify by radial position and volume
            # Substrate: r ~ (radius_inner + radius_interface)/2
            r_substrate_avg = (self.radius_inner + self.radius_interface) / 2
            # Coating: r ~ (radius_interface + radius_outer)/2
            r_coating_avg = (self.radius_interface + self.radius_outer) / 2

            if np.isclose(r_avg, r_substrate_avg, rtol=0.3):
                # Large volume near substrate radius
                if mass > 0.5 * np.pi * (self.radius_interface**2 - self.radius_inner**2) * self.height:
                    vol_substrate = tag
            elif np.isclose(r_avg, r_coating_avg, rtol=0.3):
                # Large volume near coating radius
                if mass > 0.5 * np.pi * (self.radius_outer**2 - self.radius_interface**2) * self.height:
                    vol_coating = tag
            elif np.isclose(r_avg, self.radius_interface, rtol=0.1):
                # Small volume at interface - likely disbond
                vol_disbond = tag

        # Fallback: Use ordering if classification fails
        if vol_substrate is None or vol_coating is None:
            print("Warning: Volume classification uncertain, using default ordering")
            if len(vol_tags) >= 3:
                vol_substrate = vol_tags[0]
                vol_coating = vol_tags[1]
                vol_disbond = vol_tags[2]
            else:
                raise ValueError("Fragment operation produced fewer volumes than expected")

        return (vol_substrate, vol_coating, vol_disbond)
