"""Example 02: Bonded rectangle with disbond region.

This example creates a 2D rectangle with two bonded materials and adds
a circular disbond region at the interface between them.

The disbond represents a region where the bond has failed, which will be
important for future contact mechanics and ROM analysis.
"""

from bonded_substructures.geometry import BondedRectangle
from bonded_substructures.materials.properties import ALUMINUM_7075_T6, CARBON_EPOXY_UD
from bonded_substructures.visualization import print_mesh_info


def main():
    """Create and save a bonded rectangle mesh with disbond."""

    # Define geometry parameters
    width = 10.0  # meters
    height_substrate = 2.0  # meters
    height_composite = 1.0  # meters
    depth = 2.0  # meters (thickness/depth of plate)
    mesh_size = 0.5  # meters

    # Disbond parameters (3D coordinates: x, y, z)
    disbond_center = (width / 2, height_substrate, depth / 2)  # Center of plate, at interface
    disbond_radius = 1.5  # meters

    print("Creating bonded rectangle geometry with disbond...")
    print(f"  Width: {width} m")
    print(f"  Substrate height: {height_substrate} m")
    print(f"  Composite height: {height_composite} m")
    print(f"  Mesh size: {mesh_size} m")
    print(f"  Disbond center: {disbond_center}")
    print(f"  Disbond radius: {disbond_radius} m")
    print()

    # Create geometry
    geom = BondedRectangle(
        width=width,
        height_1=height_substrate,
        height_2=height_composite,
        material_1=ALUMINUM_7075_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=mesh_size,
        depth=depth,
    )

    # Add disbond region
    print("Adding disbond region...")
    geom.add_disbond(
        position=disbond_center,
        size=disbond_radius,
        shape="circular",
    )

    # Generate mesh
    print("Generating mesh...")
    geom.generate_mesh()

    # Print mesh information
    print_mesh_info()

    # Get mesh statistics
    mesh_info = geom.get_mesh_info()
    print(f"\nMesh statistics:")
    print(f"  Nodes: {mesh_info['num_nodes']}")
    print(f"  Elements: {mesh_info['num_elements']}")
    print(f"  Has disbond: {mesh_info['has_disbond']}")

    # Save mesh
    output_file = "rectangle_disbond.msh"
    geom.save_mesh(output_file)

    # Clean up
    geom.finalize()

    print("\nDone!")
    print(f"Mesh saved to {output_file}")
    print("Use example 03 to visualize this mesh.")


if __name__ == "__main__":
    main()
