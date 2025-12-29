"""Pytest configuration and fixtures."""

import tempfile
from pathlib import Path

import pytest

from bonded_substructures.geometry import BondedRectangle, BondedCylinder
from bonded_substructures.materials.properties import ALUMINUM_6061_T6, CARBON_EPOXY_UD, ALUMINUM_7075_T6


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)


@pytest.fixture
def simple_rectangle():
    """Create a simple bonded rectangle geometry."""
    geom = BondedRectangle(
        width=10.0,
        height_1=2.0,
        height_2=1.0,
        material_1=ALUMINUM_6061_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.5,
        depth=1.0,
    )
    yield geom
    geom.finalize()


@pytest.fixture
def rectangle_with_disbond():
    """Create a bonded rectangle with a disbond region."""
    depth = 1.0
    geom = BondedRectangle(
        width=10.0,
        height_1=2.0,
        height_2=1.0,
        material_1=ALUMINUM_6061_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.5,
        depth=depth,
    )
    # 3D disbond position: (x_center, y_interface, z_center)
    geom.add_disbond(position=(5.0, 2.0, depth / 2), size=1.0, shape="circular")
    yield geom
    geom.finalize()


@pytest.fixture
def material_aluminum():
    """Aluminum material properties."""
    return ALUMINUM_6061_T6


@pytest.fixture
def material_composite():
    """Composite material properties."""
    return CARBON_EPOXY_UD


@pytest.fixture
def simple_cylinder():
    """Create a simple hollow cylinder geometry."""
    geom = BondedCylinder(
        radius_inner=0.1,
        t1=0.01,
        t2=0.02,
        height=0.5,
        material_1=ALUMINUM_7075_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.05,
    )
    yield geom
    geom.finalize()


@pytest.fixture
def cylinder_with_disbond():
    """Create hollow cylinder with disbond patch."""
    geom = BondedCylinder(
        radius_inner=0.1,
        t1=0.01,
        t2=0.02,
        height=0.5,
        material_1=ALUMINUM_7075_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.05,
    )
    # Disbond at mid-height, full circumferential band
    geom.add_disbond(
        position=(90, 0.25, 0),  # theta=90deg, z=mid-height
        size=(45, 0.1),          # 45deg arc, 0.1m axial extent
        shape="rectangular"
    )
    yield geom
    geom.finalize()
