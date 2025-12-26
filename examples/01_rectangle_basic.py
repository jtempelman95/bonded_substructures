"""Example 01: Basic bonded rectangle mesh.

This example creates a simple 2D rectangle with two bonded materials:
- Bottom half: Aluminum substrate
- Top half: Carbon/epoxy composite

The mesh is generated and saved to a file.
"""

from bonded_substructures.geometry import BondedRectangle
from bonded_substructures.materials.properties import ALUMINUM_6061_T6, CARBON_EPOXY_UD
from bonded_substructures.visualization import print_mesh_info


def main():
    """Create and save a basic bonded rectangle mesh."""

    # Define geometry parameters
    width = 10.0  # meters
    height_substrate = 2.0  # meters
    height_composite = 1.0  # meters
    depth = 2.0  # meters (thickness/depth of plate)
    mesh_size = 0.5  # meters

    print("Creating bonded 3D plate geometry...")
    print(f"  Width: {width} m")
    print(f"  Substrate height: {height_substrate} m")
    print(f"  Composite height: {height_composite} m")
    print(f"  Depth (thickness): {depth} m")
    print(f"  Mesh size: {mesh_size} m")
    print()

    # Create geometry with context manager (auto-cleanup)
    with BondedRectangle(
        width=width,
        height_1=height_substrate,
        height_2=height_composite,
        material_1=ALUMINUM_6061_T6,
        material_2=CARBON_EPOXY_UD,
        mesh_size=mesh_size,
        depth=depth,
    ) as geom:
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

        # Save mesh
        output_file = "rectangle_basic.msh"
        geom.save_mesh(output_file)

        # Optional: visualize with gmsh GUI
        # geom.visualize()

    print("\nDone!")


if __name__ == "__main__":
    main()
