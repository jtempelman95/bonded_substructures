"""Example 06: Visualize cylinder meshes with PyVista and Matplotlib.

This example demonstrates visualization of hollow cylinder meshes
created in examples 04 and 05.
"""

from pathlib import Path
import sys

# Check which visualization example to run
def main():
    """Visualize cylinder mesh with PyVista."""

    print("=" * 70)
    print("Cylinder Mesh Visualization")
    print("=" * 70)

    # Try to visualize both basic and disbond cylinders
    meshes_to_viz = [
        ("bonded_cylinder_basic.msh", "Basic Cylinder"),
        ("bonded_cylinder_disbond.msh", "Cylinder with Disbond")
    ]

    try:
        from bonded_substructures.visualization import visualize_mesh_pyvista
    except ImportError as e:
        print(f"\n⚠️  PyVista not available: {e}")
        print("Install with: conda install -c conda-forge pyvista")
        return

    for mesh_file, title in meshes_to_viz:
        if not Path(mesh_file).exists():
            print(f"\n⚠️  Mesh file not found: {mesh_file}")
            print(f"   Run the corresponding example first to generate it")
            continue

        print(f"\nVisualizing: {mesh_file}")
        print(f"Title: {title}")

        try:
            # Create output filename
            output_file = f"examples/output/mesh_{Path(mesh_file).stem}_pyvista.png"

            # Visualize
            visualize_mesh_pyvista(
                mesh_file,
                title=title,
                screenshot=output_file,
                show=False  # Don't show interactive window in batch mode
            )

            print(f"  ✓ Saved visualization: {output_file}")

        except Exception as e:
            print(f"  ✗ Visualization failed: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 70)
    print("Visualization complete!")
    print("=" * 70)


if __name__ == "__main__":
    main()
