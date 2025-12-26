"""Example 04: Create non-conformal mesh for contact mechanics.

This example demonstrates how to create a mesh with duplicate nodes at
the disbond interface, suitable for Craig-Bampton substructuring and
contact mechanics.

Steps:
1. Generate conformal mesh (with shared nodes)
2. Post-process to duplicate nodes at disbond interface
3. Result: Separate node sets for each material at interface
"""

from bonded_substructures.geometry import BondedRectangle
from bonded_substructures.materials.properties import ALUMINUM_7075_T6, CARBON_EPOXY_UD
from bonded_substructures.mesh import duplicate_disbond_interface_nodes
from bonded_substructures.visualization import print_mesh_info


def main():
    """Create non-conformal disbond mesh for substructuring."""

    # Define geometry parameters
    width = 10.0
    height_substrate = 2.0
    height_composite = 1.0
    depth = 2.0
    mesh_size = 0.5

    # Disbond parameters
    disbond_center = (width / 2, height_substrate, depth / 2)
    disbond_radius = 1.5

    print("=" * 70)
    print("Creating Non-Conformal Mesh for Craig-Bampton Substructuring")
    print("=" * 70)
    print(f"\nGeometry:")
    print(f"  Width: {width} m")
    print(f"  Substrate height: {height_substrate} m")
    print(f"  Composite height: {height_composite} m")
    print(f"  Mesh size: {mesh_size} m")
    print(f"  Disbond center: {disbond_center}")
    print(f"  Disbond radius: {disbond_radius} m\n")

    # Step 1: Create conformal mesh
    print("Step 1: Generating conformal mesh...")
    print("-" * 70)
    
    geom = BondedRectangle(
        width=width,
        height_1=height_substrate,
        height_2=height_composite,
        material_1=ALUMINUM_7075_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=mesh_size,
        depth=depth,
    )

    geom.add_disbond(
        position=disbond_center,
        size=disbond_radius,
        shape="circular",
    )

    geom.generate_mesh()
    print_mesh_info()

    # Save conformal mesh
    conformal_mesh = "rectangle_disbond_conformal.msh"
    geom.save_mesh(conformal_mesh)
    print(f"\n✓ Conformal mesh saved to: {conformal_mesh}")

    geom.finalize()

    # Step 2: Create non-conformal mesh by duplicating interface nodes
    print("\n\nStep 2: Creating non-conformal mesh...")
    print("-" * 70)
    
    nonconformal_mesh = "rectangle_disbond_nonconformal.msh"
    
    stats = duplicate_disbond_interface_nodes(
        input_mesh_file=conformal_mesh,
        output_mesh_file=nonconformal_mesh,
        disbond_physical_tag=30,      # Disbond region
        material1_physical_tag=1,     # Aluminum substrate
        material2_physical_tag=2,     # Composite coating
    )

    # Summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"Conformal mesh (shared nodes):      {conformal_mesh}")
    print(f"Non-conformal mesh (duplicate nodes): {nonconformal_mesh}")
    print(f"\nNode duplication:")
    print(f"  Material 1-Disbond interface: {stats['interface1_nodes']} nodes duplicated")
    print(f"  Material 2-Disbond interface: {stats['interface2_nodes']} nodes duplicated")
    print(f"  Total new nodes: {stats['duplicated_nodes']}")
    print(f"\nThe non-conformal mesh is ready for:")
    print(f"  - Craig-Bampton dynamic substructuring")
    print(f"  - Contact mechanics at disbond interface")
    print(f"  - Nonlinear ROM with interface DOFs")
    print("=" * 70)


if __name__ == "__main__":
    main()
