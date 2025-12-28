"""Example 01: Basic bonded plate with default dimensions.

This example creates a simple bonded plate geometry with:
- Bond interface in x-y plane (horizontal)
- Materials stacked in z-direction (through-thickness)
- Default dimensions: 1 ft × 1 ft × 0.25 in (0.15 in substrate + 0.10 in coating)

COORDINATE SYSTEM:
- x, y: In-plane dimensions (width × length)
- z: Through-thickness (material stacking)
- Bond interface at z = 0.15 in (3.81 mm)
"""

from bonded_substructures.geometry import BondedRectangle
from bonded_substructures.materials import ALUMINUM_7075_T6, CARBON_EPOXY_UD
from bonded_substructures.visualization import print_mesh_info


def main():
    """Generate basic bonded plate mesh."""

    print("=" * 70)
    print("Bonded Plate - Basic Geometry (No Disbond)")
    print("=" * 70)

    # Create bonded plate with default dimensions
    # All parameters have defaults but shown here for clarity
    plate = BondedRectangle(
        width=0.3048,           # 1 ft in x-direction
        length=0.3048,          # 1 ft in y-direction
        t1=0.00381,             # 0.15 in substrate thickness
        t2=0.00254,             # 0.10 in coating thickness
        material_1=ALUMINUM_7075_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.025,        # 25 mm (about 1 inch)
    )

    print(f"\nGeometry:")
    print(f"  Plate dimensions: {plate.width*3.28084:.2f} ft × {plate.length*3.28084:.2f} ft ({plate.width:.4f} m × {plate.length:.4f} m)")
    print(f"  Substrate (Material 1): {plate.t1*39.3701:.3f} in ({plate.t1*1000:.2f} mm) - {plate.material_1.name}")
    print(f"  Coating (Material 2): {plate.t2*39.3701:.3f} in ({plate.t2*1000:.2f} mm) - {plate.material_2.name}")
    print(f"  Total thickness: {plate.total_thickness*39.3701:.3f} in ({plate.total_thickness*1000:.2f} mm)")
    print(f"  Bond interface: x-y plane at z = {plate.t1*39.3701:.3f} in")
    print(f"  Mesh size: {plate.mesh_size*39.3701:.2f} in ({plate.mesh_size*1000:.1f} mm)")

    # Generate mesh
    print("\nGenerating mesh...")
    plate.generate_mesh()

    # Print mesh info
    print("\n" + "=" * 70)
    print_mesh_info()
    print("=" * 70)

    # Save mesh
    output_file = "bonded_plate_basic.msh"
    plate.save_mesh(output_file)
    plate.finalize()

    print(f"\n✅ Mesh saved: {output_file}")
    print("\nNext steps:")
    print("  - Run Example 02 for plate with disbond")
    print("  - Run Example 03 for visualization")


if __name__ == "__main__":
    main()
