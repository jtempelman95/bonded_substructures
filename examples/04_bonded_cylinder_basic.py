"""Example 04: Basic hollow cylinder (COPV-style) without disbond.

This example demonstrates creating a hollow cylinder geometry with radial
bonding interface, typical of Composite Overwrapped Pressure Vessels (COPVs).

The cylinder consists of:
- Inner aluminum liner (substrate)
- Outer composite overwrap (coating)
- Radial bonding interface at the cylindrical surface
"""

from pathlib import Path
from bonded_substructures.geometry import BondedCylinder
from bonded_substructures.materials import ALUMINUM_7075_T6, CARBON_EPOXY_UD
from bonded_substructures.visualization import print_mesh_info


def main():
    """Generate basic hollow cylinder mesh."""

    print("=" * 70)
    print("Hollow Cylinder - COPV Configuration (No Disbond)")
    print("=" * 70)

    # Create cylinder geometry with COPV-scale dimensions
    cylinder = BondedCylinder(
        radius_inner=0.45,          # 1.5 ft inner bore radius
        t1=0.003,                   # 3 mm aluminum liner thickness
        t2=0.020,                   # 20 mm composite overwrap thickness
        height=1.8,                 # 6 ft axial height
        material_1=ALUMINUM_7075_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=0.05,             # 50 mm characteristic element size
    )

    # Print geometry information
    print(f"\nGeometry:")
    print(f"  Coordinate system: Cylinder axis along z-direction")
    print(f"  Inner radius: {cylinder.radius_inner*3.28084:.2f} ft ({cylinder.radius_inner*1000:.1f} mm)")
    print(f"  Substrate (liner): {cylinder.material_1.name}")
    print(f"    Thickness (radial): {cylinder.t1*39.3701:.2f} in ({cylinder.t1*1000:.1f} mm)")
    print(f"    Radial extent: {cylinder.radius_inner*1000:.1f} - {cylinder.radius_interface*1000:.1f} mm")
    print(f"  Coating (overwrap): {cylinder.material_2.name}")
    print(f"    Thickness (radial): {cylinder.t2*39.3701:.2f} in ({cylinder.t2*1000:.1f} mm)")
    print(f"    Radial extent: {cylinder.radius_interface*1000:.1f} - {cylinder.radius_outer*1000:.1f} mm")
    print(f"  Interface radius: {cylinder.radius_interface*3.28084:.2f} ft ({cylinder.radius_interface*1000:.1f} mm)")
    print(f"  Outer radius: {cylinder.radius_outer*3.28084:.2f} ft ({cylinder.radius_outer*1000:.1f} mm)")
    print(f"  Height: {cylinder.height*3.28084:.2f} ft ({cylinder.height:.2f} m)")

    # Volume calculations
    import numpy as np
    vol_substrate = np.pi * (cylinder.radius_interface**2 - cylinder.radius_inner**2) * cylinder.height
    vol_coating = np.pi * (cylinder.radius_outer**2 - cylinder.radius_interface**2) * cylinder.height
    vol_total = vol_substrate + vol_coating

    print(f"\nVolumes:")
    print(f"  Substrate: {vol_substrate*1e6:.1f} cm³")
    print(f"  Coating: {vol_coating*1e6:.1f} cm³")
    print(f"  Total: {vol_total*1e6:.1f} cm³")

    # Generate mesh
    print("\nGenerating mesh...")
    print(f"  Mesh size: {cylinder.mesh_size*1000:.0f} mm")

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
    output_file = "bonded_cylinder_basic.msh"
    cylinder.save_mesh(output_file)
    cylinder.finalize()

    print(f"\n✅ Mesh saved: {output_file}")
    print(f"\nThis COPV-style cylinder geometry is ready for:")
    print(f"  • Craig-Bampton reduced order modeling")
    print(f"  • Pressure vessel analysis")
    print(f"  • Bonding interface stress analysis")
    print(f"  • Disbond studies (see Example 05)")


if __name__ == "__main__":
    main()
