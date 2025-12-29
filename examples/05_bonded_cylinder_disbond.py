"""Example 05: Hollow cylinder with circumferential disbond patch.

This example demonstrates adding a disbond patch to the cylindrical bonding
interface of a COPV-style hollow cylinder.

The disbond is specified in cylindrical coordinates:
- Circumferential position (theta in degrees)
- Axial position (z)
- Angular extent (delta_theta in degrees)
- Axial extent (delta_z)
"""

from pathlib import Path
import numpy as np
from bonded_substructures.geometry import BondedCylinder
from bonded_substructures.materials import ALUMINUM_7075_T6, CARBON_EPOXY_UD
from bonded_substructures.visualization import print_mesh_info


def main():
    """Generate hollow cylinder with disbond patch."""

    print("=" * 70)
    print("Hollow Cylinder with Circumferential Disbond Patch")
    print("=" * 70)

    # Create cylinder geometry
    cylinder = BondedCylinder(
        radius_inner=0.45,          # 1.5 ft inner bore radius
        t1=0.003,                   # 3 mm aluminum liner thickness
        t2=0.020,                   # 20 mm composite overwrap thickness
        height=1.8,                 # 6 ft axial height
        material_1=ALUMINUM_7075_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.04,             # 40 mm (finer for disbond)
    )

    # Print geometry information
    print(f"\nGeometry:")
    print(f"  Coordinate system: Cylinder axis along z-direction")
    print(f"  Inner radius: {cylinder.radius_inner*3.28084:.2f} ft ({cylinder.radius_inner*1000:.1f} mm)")
    print(f"  Interface radius: {cylinder.radius_interface*3.28084:.2f} ft ({cylinder.radius_interface*1000:.1f} mm)")
    print(f"  Outer radius: {cylinder.radius_outer*3.28084:.2f} ft ({cylinder.radius_outer*1000:.1f} mm)")
    print(f"  Height: {cylinder.height*3.28084:.2f} ft ({cylinder.height:.2f} m)")

    # Define disbond patch in cylindrical coordinates
    theta_center_deg = 90.0      # 90 degrees (along positive y-axis)
    z_center = cylinder.height / 2  # Mid-height
    theta_extent_deg = 45.0      # 45 degree arc
    z_extent = 0.2               # 0.2 m = 200 mm axial extent

    print(f"\nDisbond patch:")
    print(f"  Shape: Rectangular (in cylindrical coordinates)")
    print(f"  Circumferential position: {theta_center_deg:.0f}° (from +x axis)")
    print(f"  Axial position: {z_center:.2f} m (mid-height)")
    print(f"  Circumferential extent: {theta_extent_deg:.0f}° arc")
    print(f"  Axial extent: {z_extent*1000:.0f} mm ({z_extent*39.3701:.1f} in)")

    # Compute disbond area
    # Arc length = r * theta (radians)
    arc_length = cylinder.radius_interface * np.deg2rad(theta_extent_deg)
    disbond_area = arc_length * z_extent  # Approximate as flat patch
    interface_area = 2 * np.pi * cylinder.radius_interface * cylinder.height

    print(f"\nDisbond area:")
    print(f"  Arc length: {arc_length*1000:.0f} mm ({arc_length*39.3701:.1f} in)")
    print(f"  Patch area: {disbond_area*1550:.1f} in²")
    print(f"  Total interface area: {interface_area*1550:.1f} in²")
    print(f"  Disbond fraction: {disbond_area/interface_area*100:.1f}%")

    # Add disbond
    print(f"\nAdding disbond to geometry...")
    try:
        cylinder.add_disbond(
            position=(theta_center_deg, z_center, 0),
            size=(theta_extent_deg, z_extent),
            shape="rectangular"
        )
        print(f"  ✓ Disbond added successfully")
    except Exception as e:
        print(f"  ✗ Failed to add disbond: {e}")
        cylinder.finalize()
        return

    # Generate mesh
    print("\nGenerating mesh...")
    print(f"  Mesh size: {cylinder.mesh_size*1000:.0f} mm")
    print(f"  Note: Disbond region will be meshed")

    try:
        cylinder.generate_mesh()
        print("  ✓ Mesh generation successful")
    except Exception as e:
        print(f"  ✗ Mesh generation failed: {e}")
        import traceback
        traceback.print_exc()
        cylinder.finalize()
        return

    # Print mesh information
    print("\n" + "=" * 70)
    print_mesh_info()
    print("=" * 70)

    # Save mesh
    output_file = "bonded_cylinder_disbond.msh"
    cylinder.save_mesh(output_file)
    cylinder.finalize()

    print(f"\n✅ Mesh saved: {output_file}")

    print(f"\nKey observations:")
    print(f"  • Disbond patch created at cylindrical interface")
    print(f"  • Patch located at theta={theta_center_deg}°, z={z_center:.2f}m")
    print(f"  • This disbond can be used for:")
    print(f"    - Contact mechanics analysis (HB-AFT)")
    print(f"    - Interface stress analysis")
    print(f"    - Delamination propagation studies")
    print(f"    - COPV failure analysis")


if __name__ == "__main__":
    main()
