"""Material properties for bonded substructures."""

from dataclasses import dataclass
from typing import Optional


@dataclass
class MaterialProperties:
    """Material properties for structural analysis.

    Attributes:
        name: Material name
        E: Young's modulus (Pa) - for isotropic materials
        nu: Poisson's ratio (dimensionless)
        rho: Density (kg/m³)
        E1: Young's modulus in fiber direction (Pa) - for orthotropic materials
        E2: Young's modulus transverse to fiber (Pa) - for orthotropic materials
        E3: Young's modulus through-thickness (Pa) - for orthotropic materials
        G12: Shear modulus in 1-2 plane (Pa)
        G13: Shear modulus in 1-3 plane (Pa)
        G23: Shear modulus in 2-3 plane (Pa)
        nu12: Major Poisson's ratio
        nu13: Poisson's ratio 1-3
        nu23: Poisson's ratio 2-3
    """

    name: str
    rho: float

    # Isotropic properties
    E: Optional[float] = None
    nu: Optional[float] = None

    # Orthotropic properties
    E1: Optional[float] = None
    E2: Optional[float] = None
    E3: Optional[float] = None
    G12: Optional[float] = None
    G13: Optional[float] = None
    G23: Optional[float] = None
    nu12: Optional[float] = None
    nu13: Optional[float] = None
    nu23: Optional[float] = None

    @property
    def is_isotropic(self) -> bool:
        """Check if material is isotropic."""
        return self.E is not None and self.nu is not None

    @property
    def is_orthotropic(self) -> bool:
        """Check if material is orthotropic."""
        return all([
            self.E1 is not None,
            self.E2 is not None,
            self.E3 is not None,
            self.G12 is not None,
            self.G13 is not None,
            self.G23 is not None,
            self.nu12 is not None,
            self.nu13 is not None,
            self.nu23 is not None,
        ])

    def __post_init__(self):
        """Validate material properties."""
        if not self.is_isotropic and not self.is_orthotropic:
            raise ValueError(
                "Material must have either isotropic (E, nu) or "
                "orthotropic (E1, E2, E3, G12, G13, G23, nu12, nu13, nu23) properties"
            )


# Predefined materials
ALUMINUM_6061_T6 = MaterialProperties(
    name="Aluminum 6061-T6",
    E=68.9e9,  # Pa
    nu=0.33,
    rho=2700,  # kg/m³
)

ALUMINUM_7075_T6 = MaterialProperties(
    name="Aluminum 7075-T6",
    E=71.7e9,  # Pa
    nu=0.33,
    rho=2810,  # kg/m³
)

STEEL_4130 = MaterialProperties(
    name="Steel 4130",
    E=205e9,  # Pa
    nu=0.29,
    rho=7850,  # kg/m³
)

# Carbon fiber composite (IM7/8552) - typical unidirectional properties
# Values are representative for a unidirectional lamina
CARBON_EPOXY_UD = MaterialProperties(
    name="Carbon/Epoxy UD (IM7/8552)",
    E1=171e9,   # Fiber direction (Pa)
    E2=9.08e9,  # Transverse direction (Pa)
    E3=9.08e9,  # Through-thickness (Pa)
    G12=5.29e9, # In-plane shear (Pa)
    G13=5.29e9, # Out-of-plane shear (Pa)
    G23=2.77e9, # Transverse shear (Pa)
    nu12=0.32,  # Major Poisson's ratio
    nu13=0.32,
    nu23=0.49,
    rho=1570,   # kg/m³
)

# Glass fiber composite - typical properties
GLASS_EPOXY_UD = MaterialProperties(
    name="Glass/Epoxy UD",
    E1=45e9,    # Fiber direction (Pa)
    E2=12e9,    # Transverse direction (Pa)
    E3=12e9,    # Through-thickness (Pa)
    G12=5.5e9,  # In-plane shear (Pa)
    G13=5.5e9,  # Out-of-plane shear (Pa)
    G23=3.9e9,  # Transverse shear (Pa)
    nu12=0.28,  # Major Poisson's ratio
    nu13=0.28,
    nu23=0.40,
    rho=2000,   # kg/m³
)
