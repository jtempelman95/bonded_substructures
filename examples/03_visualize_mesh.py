"""Example 03: Visualize meshes with matplotlib and pyvista.

This example demonstrates how to visualize generated meshes using both
matplotlib (for 2D plots) and pyvista (for interactive 3D visualization).

Saves visualizations as PNG files in the examples/output directory.

Run examples 01 or 02 first to generate mesh files.
"""

import os
from pathlib import Path

import matplotlib.pyplot as plt

from bonded_substructures.visualization import plot_mesh_matplotlib, plot_mesh_pyvista


def main():
    """Visualize mesh files and save as PNG."""

    # Create output directory
    output_dir = Path("examples/output")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Mesh files to visualize
    mesh_files = {
        "rectangle_basic.msh": "basic",
        "rectangle_disbond.msh": "disbond",
    }

    print("Visualizing meshes and saving to examples/output/")
    print("=" * 60)

    for mesh_file, name in mesh_files.items():
        if not Path(mesh_file).exists():
            print(f"\n⚠️  {mesh_file} not found. Skipping.")
            print(f"   Run example 0{list(mesh_files.keys()).index(mesh_file) + 1} first.")
            continue

        print(f"\n📊 Visualizing {mesh_file}...")

        # Method 1: Matplotlib visualization (2D)
        print("   Creating matplotlib plot...")
        try:
            fig = plot_mesh_matplotlib(
                mesh_file=mesh_file,
                show_physical_groups=True,
                figsize=(14, 8),
            )

            # Save matplotlib plot
            output_path = output_dir / f"mesh_{name}_matplotlib.png"
            fig.savefig(output_path, dpi=150, bbox_inches="tight")
            print(f"   ✅ Saved: {output_path}")
            plt.close(fig)

        except FileNotFoundError:
            print(f"   ❌ Error: {mesh_file} not found.")
            continue
        except Exception as e:
            print(f"   ⚠️  Matplotlib visualization failed: {e}")

        # Method 2: PyVista visualization (save screenshot)
        print("   Creating PyVista visualization...")
        try:
            import pyvista as pv
            import meshio

            # Set PyVista to use software rendering for better WSL compatibility
            pv.OFF_SCREEN = True

            # Read mesh with meshio and convert carefully
            mesh_data = meshio.read(mesh_file)

            # Convert to PyVista UnstructuredGrid
            cells_dict = {}
            for cell_block in mesh_data.cells:
                if cell_block.type not in cells_dict:
                    cells_dict[cell_block.type] = []
                cells_dict[cell_block.type].append(cell_block.data)

            # For now, focus on tetrahedra for 3D or triangles for 2D
            mesh = None
            if "tetra" in cells_dict:
                # Create unstructured grid from tetrahedra
                import numpy as np
                cells = np.vstack(cells_dict["tetra"])
                n_cells = len(cells)
                cell_types = np.full(n_cells, pv.CellType.TETRA, dtype=np.uint8)
                vtk_cells = []
                for cell in cells:
                    vtk_cells.extend([4] + cell.tolist())
                mesh = pv.UnstructuredGrid(vtk_cells, cell_types, mesh_data.points)
            elif "triangle" in cells_dict:
                # Create polydata from triangles
                import numpy as np
                cells = np.vstack(cells_dict["triangle"])
                vtk_cells = []
                for cell in cells:
                    vtk_cells.extend([3] + cell.tolist())
                mesh = pv.PolyData(mesh_data.points, vtk_cells)

            if mesh is not None:
                # Extract surface for better visualization of 3D volumes
                if hasattr(mesh, 'extract_surface') and "tetra" in cells_dict:
                    surface = mesh.extract_surface()
                else:
                    surface = mesh

                # Create plotter with off-screen rendering
                plotter = pv.Plotter(off_screen=True, window_size=[1200, 800])

                # Add mesh with simple coloring
                plotter.add_mesh(
                    surface,
                    show_edges=True,
                    color="lightblue",
                    edge_color="black",
                    line_width=0.8,
                    opacity=0.9,
                )

                plotter.add_axes()
                plotter.show_grid()

                # Use isometric view for 3D meshes
                plotter.camera_position = 'iso'  # Valid options: xy, xz, yz, yx, zx, zy, iso
                plotter.camera.zoom(1.3)

                # Save screenshot
                output_path = output_dir / f"mesh_{name}_pyvista.png"
                plotter.screenshot(str(output_path))
                plotter.close()
                print(f"   ✅ Saved: {output_path}")
            else:
                print(f"   ⚠️  No tetrahedral or triangle cells found in mesh")

        except Exception as e:
            print(f"   ⚠️  PyVista visualization failed: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 60)
    print("✅ Visualization complete!")
    print(f"📁 Output saved to: {output_dir.absolute()}")

    # List generated files
    output_files = list(output_dir.glob("*.png"))
    if output_files:
        print("\n📄 Generated files:")
        for f in sorted(output_files):
            print(f"   - {f.name}")


if __name__ == "__main__":
    main()
