"""Example 03: Visualize meshes with matplotlib and pyvista.

This example demonstrates how to visualize generated meshes using both
matplotlib (for 2D plots) and pyvista (for interactive 3D visualization).

Run examples 01 or 02 first to generate mesh files.
"""

import matplotlib.pyplot as plt

from bonded_substructures.visualization import plot_mesh_matplotlib, plot_mesh_pyvista


def main():
    """Visualize mesh files."""

    # Choose which mesh to visualize
    mesh_files = [
        "rectangle_basic.msh",
        "rectangle_disbond.msh",
    ]

    print("Available mesh files:")
    for i, filename in enumerate(mesh_files):
        print(f"  {i+1}. {filename}")

    choice = input("\nSelect mesh file (1 or 2, or press Enter for mesh 1): ").strip()

    if choice == "2":
        mesh_file = mesh_files[1]
    else:
        mesh_file = mesh_files[0]

    print(f"\nVisualizing {mesh_file}...")

    # Method 1: Matplotlib visualization (2D)
    print("\n1. Creating matplotlib plot...")
    try:
        fig = plot_mesh_matplotlib(
            mesh_file=mesh_file,
            show_physical_groups=True,
            figsize=(12, 6),
        )
        plt.show()
    except FileNotFoundError:
        print(f"Error: {mesh_file} not found. Run example 01 or 02 first.")
        return

    # Method 2: PyVista visualization (interactive)
    print("\n2. Launching PyVista interactive viewer...")
    print("   (Close the window to continue)")
    try:
        plotter = plot_mesh_pyvista(
            mesh_file=mesh_file,
            show_edges=True,
            show_labels=False,
        )
        if plotter is not None:
            plotter.show()
    except Exception as e:
        print(f"PyVista visualization failed: {e}")
        print("This is normal if running in a headless environment.")

    print("\nDone!")


if __name__ == "__main__":
    main()
