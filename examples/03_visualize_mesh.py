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
            physical_tags = None
            cells = None

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

                # Get physical tags from cell_data (list format from meshio)
                if "gmsh:physical" in mesh_data.cell_data:
                    # Find tetra blocks and concatenate their physical tags
                    phys_tags_blocks = []
                    for i, cell_block in enumerate(mesh_data.cells):
                        if cell_block.type == "tetra":
                            phys_tags_blocks.append(mesh_data.cell_data["gmsh:physical"][i])
                    if phys_tags_blocks:
                        physical_tags = np.concatenate(phys_tags_blocks)
                        mesh.cell_data["Material"] = physical_tags
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

                    # Transfer material data from volume to surface based on geometry
                    # Color surface elements by their y-coordinate
                    if mesh.n_points > 0:
                        import numpy as np

                        # Find actual interface height from the volume mesh physical tags
                        # Get y-coordinates of nodes in each material
                        mat1_y_coords = []
                        mat2_y_coords = []

                        if physical_tags is not None and cells is not None:
                            for i, tag in enumerate(physical_tags):
                                elem = cells[i]  # Tetrahedral element
                                elem_y_coords = [mesh.points[node, 1] for node in elem]
                                if tag == 1:  # Material 1
                                    mat1_y_coords.extend(elem_y_coords)
                                elif tag == 2:  # Material 2
                                    mat2_y_coords.extend(elem_y_coords)

                        # Interface is at top of Material 1 / bottom of Material 2
                        if mat1_y_coords and mat2_y_coords:
                            mat1_max_y = max(mat1_y_coords)
                            mat2_min_y = min(mat2_y_coords)
                            y_interface = (mat1_max_y + mat2_min_y) / 2
                        else:
                            y_coords = surface.points[:, 1]
                            y_min, y_max = y_coords.min(), y_coords.max()
                            y_interface = (y_min + y_max) / 2

                        y_coords = surface.points[:, 1]
                        y_min, y_max = y_coords.min(), y_coords.max()
                        interface_tolerance = (y_max - y_min) * 0.01  # Very thin interface (1% of height)

                        # Assign material IDs based on cell centroid y-position
                        material_ids = np.zeros(surface.n_cells)
                        for i in range(surface.n_cells):
                            cell = surface.get_cell(i)
                            point_ids = cell.point_ids
                            centroid_y = np.mean([surface.points[pid, 1] for pid in point_ids])

                            if abs(centroid_y - y_interface) < interface_tolerance:
                                material_ids[i] = 30  # Interface
                            elif centroid_y < y_interface:
                                material_ids[i] = 1   # Material 1
                            else:
                                material_ids[i] = 2   # Material 2

                        surface.cell_data["Material"] = material_ids
                else:
                    surface = mesh

                # Create plotter with off-screen rendering
                plotter = pv.Plotter(off_screen=True, window_size=[1200, 800])

                # Add mesh with material-based coloring
                if "Material" in surface.array_names:
                    # Create custom colormap for discrete material IDs
                    import matplotlib.pyplot as plt
                    from matplotlib.colors import ListedColormap

                    # Map material IDs to colors
                    colors_dict = {
                        1: [1.0, 0.42, 0.42],   # Red for Material 1
                        2: [0.31, 0.80, 0.77],  # Cyan for Material 2
                        30: [1.0, 0.90, 0.43],  # Yellow for Disbond
                    }

                    # Get unique material IDs and assign colors
                    unique_ids = np.unique(surface["Material"])
                    color_array = np.zeros((surface.n_cells, 3))
                    for i, mat_id in enumerate(surface["Material"]):
                        if mat_id in colors_dict:
                            color_array[i] = colors_dict[mat_id]

                    plotter.add_mesh(
                        surface,
                        show_edges=True,
                        rgb=True,
                        scalars=color_array,
                        edge_color="black",
                        line_width=0.5,
                        opacity=0.95,
                    )

                    # Add custom legend
                    plotter.add_text(
                        "Material 1 (Red)\nMaterial 2 (Cyan)\nDisbond (Yellow)",
                        position='upper_right',
                        font_size=10,
                        color='black'
                    )
                else:
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
