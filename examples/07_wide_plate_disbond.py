"""Example 07: Wide plate with localized disbond region.

This example creates a wider plate geometry with a small, localized disbond
for realistic structural analysis scenarios.

COORDINATE SYSTEM:
- x-y plane: In-plane dimensions (width × length)
- z-direction: Through-thickness (material stacking)
- Bond interface: x-y plane at z = t1

Geometry:
- Width (x): 20 m (wide plate)
- Length (y): 2 m
- Substrate thickness (t1): 2.0 m (aluminum)
- Coating thickness (t2): 1.0 m (carbon/epoxy)
- Disbond: Small circular region (radius 0.75m) at center of interface
- Mesh size: 0.3 m (finer mesh for better resolution)
"""

import gmsh

from bonded_substructures.geometry import BondedRectangle
from bonded_substructures.materials import ALUMINUM_7075_T6, CARBON_EPOXY_UD
from bonded_substructures.visualization import print_mesh_info


def main():
    """Generate wide plate with localized disbond."""

    print("=" * 70)
    print("Wide Plate with Localized Disbond")
    print("=" * 70)

    # Create wide plate geometry
    width = 20.0  # Wide plate in x-direction
    length = 2.0  # Length in y-direction
    t1 = 2.0  # Aluminum substrate thickness (through-thickness)
    t2 = 1.0  # Carbon/Epoxy coating thickness
    mesh_size = 0.3  # Finer mesh for better resolution

    print(f"\nGeometry Parameters:")
    print(f"  Width (x): {width} m (wide plate)")
    print(f"  Length (y): {length} m")
    print(f"  Substrate thickness (t1): {t1} m")
    print(f"  Coating thickness (t2): {t2} m")
    print(f"  Total thickness: {t1 + t2} m")
    print(f"  Bond interface: x-y plane at z = {t1} m")
    print(f"  Mesh size: {mesh_size} m")

    # Create rectangle with disbond
    rect = BondedRectangle(
        width=width,
        length=length,
        t1=t1,
        t2=t2,
        material_1=ALUMINUM_7075_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=mesh_size,
    )

    # Add small localized disbond at center
    disbond_radius = 0.75  # Small, localized disbond
    disbond_center = (width/2, length/2, t1)  # Center of x-y plane, at interface z = t1

    print(f"\nDisbond Parameters:")
    print(f"  Radius: {disbond_radius} m (localized)")
    print(f"  Center: ({disbond_center[0]}, {disbond_center[1]}, {disbond_center[2]} m)")
    print(f"  Location: x-y plane at bond interface (z = {t1} m)")
    print(f"  Disbond area: {3.14159 * disbond_radius**2:.2f} m²")
    print(f"  Total interface area: {width * length:.2f} m²")
    print(f"  Disbond fraction: {(3.14159 * disbond_radius**2) / (width * length) * 100:.1f}%")

    rect.add_disbond(
        position=disbond_center,
        size=disbond_radius,
        shape="circular"
    )

    # Generate mesh
    print("\nGenerating mesh...")
    print(f"  This may take a moment for the wide plate geometry...")

    rect.generate_mesh()

    # Print mesh info
    print("\n" + "=" * 70)
    print_mesh_info()
    print("=" * 70)

    # Save mesh
    mesh_file = "wide_plate_disbond.msh"
    rect.save_mesh(mesh_file)
    rect.finalize()

    print(f"\n✅ Wide plate mesh generated: {mesh_file}")
    print("\nNext steps:")
    print("  - Run Example 08 for time-domain response simulation")
    print("  - Use this mesh for Craig-Bampton ROM analysis")


if __name__ == "__main__":
    main()
