"""Tests for mesh markers and physical tags."""

import pytest

from bonded_substructures.mesh.markers import (
    PhysicalTag,
    get_boundary_tags,
    get_interface_tags,
    get_material_tags,
    get_tag_name,
)


def test_physical_tag_values():
    """Test that physical tag values are correct."""
    assert PhysicalTag.MATERIAL_1 == 1
    assert PhysicalTag.MATERIAL_2 == 2
    assert PhysicalTag.INTERFACE_BONDED == 10
    assert PhysicalTag.BOUNDARY_LEFT == 20


def test_get_material_tags():
    """Test getting all material tags."""
    material_tags = get_material_tags()
    assert PhysicalTag.MATERIAL_1 in material_tags
    assert PhysicalTag.MATERIAL_2 in material_tags
    assert len(material_tags) == 3


def test_get_interface_tags():
    """Test getting all interface tags."""
    interface_tags = get_interface_tags()
    assert PhysicalTag.INTERFACE_BONDED in interface_tags
    assert PhysicalTag.INTERFACE_DISBOND in interface_tags
    assert len(interface_tags) == 2


def test_get_boundary_tags():
    """Test getting all boundary tags."""
    boundary_tags = get_boundary_tags()
    assert PhysicalTag.BOUNDARY_LEFT in boundary_tags
    assert PhysicalTag.BOUNDARY_RIGHT in boundary_tags
    assert PhysicalTag.BOUNDARY_TOP in boundary_tags
    assert PhysicalTag.BOUNDARY_BOTTOM in boundary_tags


def test_get_tag_name():
    """Test getting tag name from integer."""
    assert get_tag_name(1) == "MATERIAL_1"
    assert get_tag_name(10) == "INTERFACE_BONDED"
    assert get_tag_name(20) == "BOUNDARY_LEFT"


def test_get_tag_name_unknown():
    """Test getting name for unknown tag."""
    name = get_tag_name(999)
    assert name.startswith("UNKNOWN_TAG_")
    assert "999" in name


def test_tag_enum_membership():
    """Test that tags are valid enum members."""
    assert PhysicalTag.MATERIAL_1 in PhysicalTag
    assert PhysicalTag.DISBOND_REGION in PhysicalTag
    assert 999 not in PhysicalTag
