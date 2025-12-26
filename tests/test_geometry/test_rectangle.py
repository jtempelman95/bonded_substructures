"""Tests for BondedRectangle geometry."""

import pytest
import gmsh

from bonded_substructures.geometry import BondedRectangle
from bonded_substructures.materials.properties import ALUMINUM_6061_T6, CARBON_EPOXY_UD
from bonded_substructures.mesh.markers import PhysicalTag


def test_rectangle_creation(simple_rectangle):
    """Test that a rectangle geometry can be created."""
    assert simple_rectangle is not None
    assert simple_rectangle.width == 10.0
    assert simple_rectangle.height_1 == 2.0
    assert simple_rectangle.height_2 == 1.0
    assert simple_rectangle.total_height == 3.0


def test_rectangle_materials(simple_rectangle):
    """Test that materials are properly assigned."""
    assert simple_rectangle.material_1.name == "Aluminum 6061-T6"
    assert simple_rectangle.material_2.name == "Carbon/Epoxy UD (IM7/8552)"


def test_mesh_generation(simple_rectangle):
    """Test that mesh can be generated."""
    simple_rectangle.generate_mesh()

    # Check that gmsh is initialized
    assert gmsh.isInitialized()

    # Check mesh info
    mesh_info = simple_rectangle.get_mesh_info()
    assert mesh_info["num_nodes"] > 0
    assert mesh_info["num_elements"] > 0
    assert mesh_info["has_disbond"] is False


def test_physical_groups(simple_rectangle):
    """Test that physical groups are created correctly."""
    simple_rectangle.generate_mesh()

    physical_groups = simple_rectangle.get_physical_tags()

    # Should have material and boundary groups
    assert len(physical_groups) > 0
    assert "Aluminum 6061-T6" in physical_groups
    assert "Carbon/Epoxy UD (IM7/8552)" in physical_groups


def test_disbond_addition():
    """Test adding a disbond region."""
    geom = BondedRectangle(
        width=10.0,
        height_1=2.0,
        height_2=1.0,
        material_1=ALUMINUM_6061_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.5,
    )

    # Initially no disbond
    assert geom._has_disbond is False

    # Add disbond
    geom.add_disbond(position=(5.0, 2.0), size=1.0, shape="circular")

    assert geom._has_disbond is True
    assert geom._disbond_position == (5.0, 2.0)
    assert geom._disbond_size == 1.0
    assert geom._disbond_shape == "circular"

    geom.finalize()


def test_disbond_mesh_generation(rectangle_with_disbond):
    """Test mesh generation with disbond."""
    rectangle_with_disbond.generate_mesh()

    mesh_info = rectangle_with_disbond.get_mesh_info()
    assert mesh_info["has_disbond"] is True
    assert mesh_info["num_nodes"] > 0
    assert mesh_info["num_elements"] > 0


def test_invalid_disbond_shape():
    """Test that invalid disbond shape raises error."""
    geom = BondedRectangle(
        width=10.0,
        height_1=2.0,
        height_2=1.0,
        material_1=ALUMINUM_6061_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.5,
    )

    with pytest.raises(ValueError):
        geom.add_disbond(position=(5.0, 2.0), size=1.0, shape="invalid")

    geom.finalize()


def test_mesh_save(simple_rectangle, temp_dir):
    """Test saving mesh to file."""
    simple_rectangle.generate_mesh()

    output_file = temp_dir / "test_mesh.msh"
    simple_rectangle.save_mesh(str(output_file))

    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_context_manager():
    """Test using geometry as context manager."""
    with BondedRectangle(
        width=10.0,
        height_1=2.0,
        height_2=1.0,
        material_1=ALUMINUM_6061_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.5,
    ) as geom:
        geom.generate_mesh()
        assert gmsh.isInitialized()

    # After exiting context, gmsh should be finalized
    # (Note: gmsh.isInitialized() may still return True if another test initialized it)
