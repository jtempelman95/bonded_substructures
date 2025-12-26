"""Physical tag definitions for mesh regions."""

from enum import IntEnum


class PhysicalTag(IntEnum):
    """Physical tags for marking mesh regions.

    These tags are used by gmsh to identify different regions of the mesh.
    Tags are organized by category:
    - 1-9: Material regions
    - 10-19: Interfaces
    - 20-29: Boundaries
    - 30-39: Special regions (disbonds, etc.)
    """

    # Material regions
    MATERIAL_1 = 1  # Primary material (e.g., substrate, aluminum)
    MATERIAL_2 = 2  # Secondary material (e.g., coating, composite)
    MATERIAL_3 = 3  # Tertiary material (for multi-layer structures)

    # Interface regions
    INTERFACE_BONDED = 10      # Bonded interface between materials
    INTERFACE_DISBOND = 11     # Disbond interface (contact region)

    # Boundary markers (2D rectangle)
    BOUNDARY_LEFT = 20
    BOUNDARY_RIGHT = 21
    BOUNDARY_BOTTOM = 22
    BOUNDARY_TOP = 23

    # Boundary markers (3D cylinder) - for future use
    BOUNDARY_INNER = 24  # Inner surface
    BOUNDARY_OUTER = 25  # Outer surface
    BOUNDARY_END_1 = 26  # End cap 1
    BOUNDARY_END_2 = 27  # End cap 2

    # Special regions
    DISBOND_REGION = 30  # Volume/area of disbond


def get_material_tags():
    """Get all material region tags."""
    return [
        PhysicalTag.MATERIAL_1,
        PhysicalTag.MATERIAL_2,
        PhysicalTag.MATERIAL_3,
    ]


def get_interface_tags():
    """Get all interface tags."""
    return [
        PhysicalTag.INTERFACE_BONDED,
        PhysicalTag.INTERFACE_DISBOND,
    ]


def get_boundary_tags():
    """Get all boundary tags."""
    return [
        PhysicalTag.BOUNDARY_LEFT,
        PhysicalTag.BOUNDARY_RIGHT,
        PhysicalTag.BOUNDARY_BOTTOM,
        PhysicalTag.BOUNDARY_TOP,
        PhysicalTag.BOUNDARY_INNER,
        PhysicalTag.BOUNDARY_OUTER,
        PhysicalTag.BOUNDARY_END_1,
        PhysicalTag.BOUNDARY_END_2,
    ]


def get_tag_name(tag: int) -> str:
    """Get human-readable name for a tag.

    Args:
        tag: Physical tag integer value

    Returns:
        Name of the tag
    """
    try:
        return PhysicalTag(tag).name
    except ValueError:
        return f"UNKNOWN_TAG_{tag}"
