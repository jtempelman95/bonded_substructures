"""Example 02: Bonded plate with localized disbond.

This example creates a bonded plate with a small circular disbond
at the bond interface (x-y plane at z = t1).

Default geometry:
- Plate: 1 ft × 1 ft × 0.25 in
- Disbond: 1 in radius circular region at center
- Interface: x-y plane at z = 0.15 in (3.81 mm)
"""

from bonded_substructures.geometry import BondedRectangle
from bonded_substructures.materials import ALUMINUM_7075_T6, CARBON_EPOXY_UD
from bonded_substructures.visualization import print_mesh_info


def main():
    """Generate bonded plate with disbond."""

    print("=" * 70)
    print("Bonded Plate with Localized Disbond")
    print("=" * 70)

    # Create bonded plate with default dimensions
    plate = BondedRectangle(
        width=0.3048,           # 1 ft
        length=0.3048,          # 1 ft
        t1=0.00381,             # 0.15 in substrate
        t2=0.00254,             # 0.10 in coating
        material_1=ALUMINUM_7075_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.0127,       # 0.5 in mesh (12.7 mm)
    )

    print(f"\nPlate geometry:")
    print(f"  Dimensions: {plate.width*3.28084:.2f} ft × {plate.length*3.28084:.2f} ft ({plate.width:.4f} m × {plate.length:.4f} m)")
    print(f"  Substrate: {plate.t1*39.3701:.3f} in ({plate.t1*1000:.2f} mm) - {plate.material_1.name}")
    print(f"  Coating: {plate.t2*39.3701:.3f} in ({plate.t2*1000:.2f} mm) - {plate.material_2.name}")
    print(f"  Total thickness: {plate.total_thickness*39.3701:.3f} in ({plate.total_thickness*1000:.2f} mm)")

    # Add circular disbond at center of plate
    disbond_radius = 0.0254  # 1 in radius
    disbond_center = (plate.width/2, plate.length/2, plate.t1)  # Center, at interface

    print(f"\nDisbond:")
    print(f"  Shape: Circular")
    print(f"  Radius: {disbond_radius*39.3701:.2f} in ({disbond_radius*1000:.1f} mm)")
    print(f"  Center: ({disbond_center[0]*3.28084:.2f} ft, {disbond_center[1]*3.28084:.2f} ft, {disbond_center[2]*39.3701:.3f} in)")
    print(f"  Location: x-y plane at bond interface (z = {plate.t1*39.3701:.3f} in)")

    # Calculate disbond area fraction
    disbond_area = 3.14159 * disbond_radius**2
    total_area = plate.width * plate.length
    print(f"  Disbond area: {disbond_area*1550:.2f} in² ({disbond_area/total_area*100:.1f}% of interface)")

    plate.add_disbond(
        position=disbond_center,
        size=disbond_radius,
        shape="circular"
    )

    # Generate mesh
    print("\nGenerating mesh...")
    plate.generate_mesh()

    # Print mesh info
    print("\n" + "=" * 70)
    print_mesh_info()
    print("=" * 70)

    # Save mesh
    output_file = "bonded_plate_disbond.msh"
    plate.save_mesh(output_file)
    plate.finalize()

    print(f"\n✅ Mesh saved: {output_file}")
    print("\nNext steps:")
    print("  - Run Example 03 for visualization")
    print("  - Run Example 04 for Craig-Bampton ROM")
    print("  - Run Example 05 for time-domain response")


if __name__ == "__main__":
    main()
