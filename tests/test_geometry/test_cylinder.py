"""Tests for BondedCylinder geometry."""

import pytest
import numpy as np

from bonded_substructures.geometry import BondedCylinder
from bonded_substructures.materials import ALUMINUM_7075_T6, CARBON_EPOXY_UD


def test_cylinder_creation(simple_cylinder):
    """Test that a cylinder geometry can be created."""
    assert simple_cylinder is not None
    assert simple_cylinder.radius_inner == 0.1
    assert simple_cylinder.t1 == 0.01
    assert simple_cylinder.t2 == 0.02
    assert simple_cylinder.height == 0.5


def test_cylinder_computed_radii(simple_cylinder):
    """Test that computed radii are correct."""
    assert np.isclose(simple_cylinder.radius_interface, 0.11)  # 0.1 + 0.01
    assert np.isclose(simple_cylinder.radius_outer, 0.13)      # 0.1 + 0.01 + 0.02


def test_cylinder_materials(simple_cylinder):
    """Test that materials are properly assigned."""
    assert simple_cylinder.material_1.name == "Aluminum 7075-T6"
    assert simple_cylinder.material_2.name == "Carbon/Epoxy UD (IM7/8552)"


def test_cylinder_mesh_generation(simple_cylinder):
    """Test that mesh can be generated."""
    simple_cylinder.generate_mesh()
    # If we get here without error, mesh generation succeeded
    assert True


def test_cylinder_disbond_addition():
    """Test adding a disbond region."""
    geom = BondedCylinder(
        radius_inner=0.1,
        t1=0.01,
        t2=0.02,
        height=0.5,
        mesh_size=0.05,
    )

    assert geom._has_disbond is False

    geom.add_disbond(
        position=(90, 0.25, 0),
        size=(45, 0.1),
        shape="rectangular"
    )

    assert geom._has_disbond is True
    assert geom._disbond_position == (90, 0.25, 0)
    assert geom._disbond_size == (45, 0.1)
    assert geom._disbond_shape == "rectangular"

    geom.finalize()


def test_cylinder_invalid_disbond_shape():
    """Test that invalid disbond shape raises error."""
    geom = BondedCylinder(
        radius_inner=0.1,
        t1=0.01,
        t2=0.02,
        height=0.5,
        mesh_size=0.05,
    )

    with pytest.raises(ValueError, match="only supports 'rectangular'"):
        geom.add_disbond(
            position=(90, 0.25, 0),
            size=(45, 0.1),
            shape="circular"
        )

    geom.finalize()


def test_cylinder_disbond_validation():
    """Test disbond position validation."""
    geom = BondedCylinder(
        radius_inner=0.1,
        t1=0.01,
        t2=0.02,
        height=0.5,
        mesh_size=0.05,
    )

    # z-position outside cylinder height should fail
    with pytest.raises(ValueError, match="outside cylinder height"):
        geom.add_disbond(
            position=(90, 1.0, 0),  # z=1.0 > height=0.5
            size=(45, 0.1),
            shape="rectangular"
        )

    geom.finalize()


def test_cylinder_context_manager():
    """Test using cylinder as context manager."""
    with BondedCylinder(
        radius_inner=0.1,
        t1=0.01,
        t2=0.02,
        height=0.5,
        mesh_size=0.05,
    ) as geom:
        geom.generate_mesh()
        assert geom is not None

    # After exiting context, gmsh should be finalized
    # No assertion here - just verify no errors
