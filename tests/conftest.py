"""Pytest configuration and fixtures."""

import tempfile
from pathlib import Path

import pytest

from bonded_substructures.geometry import BondedRectangle
from bonded_substructures.materials.properties import ALUMINUM_6061_T6, CARBON_EPOXY_UD


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
