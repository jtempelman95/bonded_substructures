"""Mesh utilities and physical tag management."""

from bonded_substructures.mesh.duplicate_nodes import duplicate_disbond_interface_nodes
from bonded_substructures.mesh.markers import PhysicalTag

__all__ = ["PhysicalTag", "duplicate_disbond_interface_nodes"]
