"""Create PyVista visualizations of COPV with disbond region highlighted."""

import numpy as np
from pathlib import Path

try:
    import pyvista as pv
    import meshio
    PYVISTA_AVAILABLE = True
except ImportError as e:
    print(f"Error: PyVista or meshio not available: {e}")
    print("Install with: conda install -c conda-forge pyvista meshio")
    exit(1)


def load_gmsh_mesh_to_pyvista(mesh_file):
    """Load gmsh mesh and convert to PyVista with proper cell data.

    Args:
        mesh_file: Path to .msh file

    Returns:
        pv.UnstructuredGrid with material tags
    """
    # Read mesh with meshio
    mesh = meshio.read(mesh_file)

    # Extract points
    points = mesh.points

    # Collect ALL tetrahedral cells and their data from all blocks
    all_tetra_cells = []
    all_tetra_data = []

    for i, cell_block in enumerate(mesh.cells):
        if cell_block.type == "tetra":
            all_tetra_cells.append(cell_block.data)
            # Get corresponding physical tags
            if "gmsh:physical" in mesh.cell_data:
                all_tetra_data.append(mesh.cell_data["gmsh:physical"][i])

    if not all_tetra_cells:
        raise ValueError("No tetrahedral cells found in mesh")

    # Concatenate all tetrahedral cells
    tetra_cells = np.vstack(all_tetra_cells)

    # Concatenate all material data
    if all_tetra_data:
        tetra_data = np.concatenate(all_tetra_data)
    else:
        tetra_data = None

    # Create PyVista UnstructuredGrid
    # Format: [4, p0, p1, p2, p3, 4, p0, p1, p2, p3, ...]
    cells = np.hstack([np.full((tetra_cells.shape[0], 1), 4), tetra_cells]).ravel()
    cell_types = np.full(tetra_cells.shape[0], pv.CellType.TETRA)

    grid = pv.UnstructuredGrid(cells, cell_types, points)

    # Add material data if available
    if tetra_data is not None:
        grid["Material"] = tetra_data

    return grid


def create_copv_visualization(mesh_file, title, output_dir, config_name):
    """Create comprehensive PyVista visualization of COPV.

    Args:
        mesh_file: Path to mesh file
        title: Plot title
        output_dir: Directory for output images
        config_name: Name for output files (e.g., 'ideal', 'disbond')
    """
    print(f"\nProcessing: {title}")
    print(f"Mesh file: {mesh_file}")

    # Load mesh
    grid = load_gmsh_mesh_to_pyvista(mesh_file)

    print(f"  Loaded: {grid.n_points:,} points, {grid.n_cells:,} cells")

    # Check for material tags
    has_materials = "Material" in grid.array_names
    if has_materials:
        unique_materials = np.unique(grid["Material"])
        print(f"  Material tags: {unique_materials}")
        has_disbond = 30 in unique_materials  # Disbond tag is 30
    else:
        print("  Warning: No material tags found")
        has_disbond = False

    # Create visualizations
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. Full mesh view with materials
    print("  Creating full mesh view...")
    pl = pv.Plotter(off_screen=True, window_size=[1920, 1080])

    if has_materials:
        pl.add_mesh(
            grid,
            scalars="Material",
            show_edges=False,
            cmap="Set1",
            show_scalar_bar=True,
            scalar_bar_args={
                'title': 'Material Tag',
                'title_font_size': 20,
                'label_font_size': 16,
                'n_labels': 5,
                'position_x': 0.85,
                'position_y': 0.3,
            }
        )
    else:
        pl.add_mesh(grid, color='lightblue', show_edges=True)

    pl.camera_position = 'iso'
    pl.add_title(title, font_size=20)
    pl.add_axes(
        xlabel='X (m)',
        ylabel='Y (m)',
        zlabel='Z (m)',
        line_width=5,
        labels_off=False
    )
    pl.background_color = 'white'

    output_file = output_dir / f"copv_{config_name}_full_mesh.png"
    pl.screenshot(str(output_file), return_img=False)
    pl.close()
    print(f"    ✓ Saved: {output_file}")

    # 2. Disbond highlighted view (if disbond present)
    if has_disbond:
        print("  Creating disbond-highlighted view...")

        # Extract disbond region
        disbond_mask = grid["Material"] == 30
        disbond_grid = grid.extract_cells(disbond_mask)

        # Extract substrate and coating
        substrate_mask = grid["Material"] == 1
        coating_mask = grid["Material"] == 2
        substrate_grid = grid.extract_cells(substrate_mask)
        coating_grid = grid.extract_cells(coating_mask)

        print(f"    Disbond cells: {disbond_grid.n_cells:,}")
        print(f"    Substrate cells: {substrate_grid.n_cells:,}")
        print(f"    Coating cells: {coating_grid.n_cells:,}")

        # Create plot with disbond highlighted
        pl = pv.Plotter(off_screen=True, window_size=[1920, 1080])

        # Add substrate (semi-transparent)
        pl.add_mesh(
            substrate_grid,
            color='steelblue',
            opacity=0.3,
            show_edges=False,
            label='Substrate (Al 7075-T6)'
        )

        # Add coating (semi-transparent)
        pl.add_mesh(
            coating_grid,
            color='lightgreen',
            opacity=0.3,
            show_edges=False,
            label='Coating (C/E)'
        )

        # Add disbond (opaque, bright color)
        pl.add_mesh(
            disbond_grid,
            color='red',
            opacity=1.0,
            show_edges=True,
            edge_color='darkred',
            line_width=2,
            label='Disbond Region'
        )

        pl.add_legend(size=(0.3, 0.3), loc='upper right')
        pl.camera_position = 'iso'
        pl.add_title(f"{title}\nDisbond Region Highlighted", font_size=20)
        pl.add_axes(xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)')
        pl.background_color = 'white'

        output_file = output_dir / f"copv_{config_name}_disbond_highlighted.png"
        pl.screenshot(str(output_file), return_img=False)
        pl.close()
        print(f"    ✓ Saved: {output_file}")

        # 3. Disbond only (closeup)
        print("  Creating disbond closeup...")
        pl = pv.Plotter(off_screen=True, window_size=[1920, 1080])

        pl.add_mesh(
            disbond_grid,
            color='red',
            show_edges=True,
            edge_color='black',
            line_width=1,
        )

        pl.camera_position = 'iso'
        pl.add_title(f"Disbond Region Only\nCircumferential Band at Interface", font_size=20)
        pl.add_axes(xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)')
        pl.background_color = 'white'

        output_file = output_dir / f"copv_{config_name}_disbond_only.png"
        pl.screenshot(str(output_file), return_img=False)
        pl.close()
        print(f"    ✓ Saved: {output_file}")

        # 4. Cross-section through disbond
        print("  Creating cross-section view...")

        # Create slice through middle of disbond
        disbond_center = disbond_grid.center

        # Slice in XZ plane (through y=0)
        slice_xz = grid.slice(normal='y', origin=disbond_center)

        pl = pv.Plotter(off_screen=True, window_size=[1920, 1080])

        if "Material" in slice_xz.array_names:
            pl.add_mesh(
                slice_xz,
                scalars="Material",
                cmap="Set1",
                show_edges=True,
                edge_color='black',
                line_width=1,
                show_scalar_bar=True,
                scalar_bar_args={'title': 'Material Tag'}
            )
        else:
            pl.add_mesh(slice_xz, color='lightblue', show_edges=True)

        # Add outline
        pl.add_mesh(grid.outline(), color='black', line_width=3)

        pl.camera_position = 'xz'
        pl.add_title(f"{title}\nCross-Section (XZ Plane)", font_size=20)
        pl.add_axes(xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)')
        pl.background_color = 'white'

        output_file = output_dir / f"copv_{config_name}_cross_section.png"
        pl.screenshot(str(output_file), return_img=False)
        pl.close()
        print(f"    ✓ Saved: {output_file}")

    else:
        # For ideal COPV, show material separation
        print("  Creating material separation view...")

        substrate_mask = grid["Material"] == 1
        coating_mask = grid["Material"] == 2
        substrate_grid = grid.extract_cells(substrate_mask)
        coating_grid = grid.extract_cells(coating_mask)

        pl = pv.Plotter(off_screen=True, window_size=[1920, 1080])

        pl.add_mesh(
            substrate_grid,
            color='steelblue',
            opacity=0.7,
            show_edges=False,
            label='Substrate (Al 7075-T6)'
        )

        pl.add_mesh(
            coating_grid,
            color='lightgreen',
            opacity=0.7,
            show_edges=False,
            label='Coating (C/E)'
        )

        pl.add_legend(size=(0.3, 0.2), loc='upper right')
        pl.camera_position = 'iso'
        pl.add_title(f"{title}\nMaterial Separation", font_size=20)
        pl.add_axes(xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)')
        pl.background_color = 'white'

        output_file = output_dir / f"copv_{config_name}_materials.png"
        pl.screenshot(str(output_file), return_img=False)
        pl.close()
        print(f"    ✓ Saved: {output_file}")

    # 5. Multiple viewpoints montage
    print("  Creating multi-view montage...")

    pl = pv.Plotter(shape=(2, 2), off_screen=True, window_size=[1920, 1920])

    # View 1: Isometric
    pl.subplot(0, 0)
    if has_materials:
        pl.add_mesh(grid, scalars="Material", cmap="Set1", show_edges=False)
    else:
        pl.add_mesh(grid, color='lightblue')
    pl.camera_position = 'iso'
    pl.add_title("Isometric View", font_size=16)
    pl.add_axes(line_width=3)

    # View 2: Top (XY)
    pl.subplot(0, 1)
    if has_materials:
        pl.add_mesh(grid, scalars="Material", cmap="Set1", show_edges=False)
    else:
        pl.add_mesh(grid, color='lightblue')
    pl.camera_position = 'xy'
    pl.add_title("Top View (XY)", font_size=16)
    pl.add_axes(line_width=3)

    # View 3: Side (XZ)
    pl.subplot(1, 0)
    if has_materials:
        pl.add_mesh(grid, scalars="Material", cmap="Set1", show_edges=False)
    else:
        pl.add_mesh(grid, color='lightblue')
    pl.camera_position = 'xz'
    pl.add_title("Side View (XZ)", font_size=16)
    pl.add_axes(line_width=3)

    # View 4: Front (YZ)
    pl.subplot(1, 1)
    if has_materials:
        pl.add_mesh(grid, scalars="Material", cmap="Set1", show_edges=False)
    else:
        pl.add_mesh(grid, color='lightblue')
    pl.camera_position = 'yz'
    pl.add_title("Front View (YZ)", font_size=16)
    pl.add_axes(line_width=3)

    pl.background_color = 'white'

    output_file = output_dir / f"copv_{config_name}_multiview.png"
    pl.screenshot(str(output_file), return_img=False)
    pl.close()
    print(f"    ✓ Saved: {output_file}")


def main():
    """Create all PyVista COPV visualizations."""

    print("=" * 70)
    print("COPV PyVista Visualization Suite")
    print("=" * 70)

    output_dir = Path("examples/output")

    # Visualize both configurations
    configs = [
        {
            "file": "bonded_cylinder_basic.msh",
            "title": "Ideal COPV (No Disbond)",
            "name": "ideal"
        },
        {
            "file": "bonded_cylinder_disbond.msh",
            "title": "COPV with Circumferential Disbond",
            "name": "disbond"
        }
    ]

    for config in configs:
        mesh_file = Path(config["file"])

        if not mesh_file.exists():
            print(f"\n⚠️  Mesh file not found: {mesh_file}")
            print(f"   Run examples to generate meshes first")
            continue

        try:
            create_copv_visualization(
                mesh_file,
                config["title"],
                output_dir,
                config["name"]
            )
        except Exception as e:
            print(f"\n✗ Error processing {mesh_file}: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 70)
    print("PyVista Visualization Complete!")
    print("=" * 70)
    print(f"\nGenerated files in: {output_dir}/")
    print("\nIdeal COPV:")
    print("  - copv_ideal_full_mesh.png")
    print("  - copv_ideal_materials.png")
    print("  - copv_ideal_multiview.png")
    print("\nDisbonded COPV:")
    print("  - copv_disbond_full_mesh.png")
    print("  - copv_disbond_disbond_highlighted.png")
    print("  - copv_disbond_disbond_only.png")
    print("  - copv_disbond_cross_section.png")
    print("  - copv_disbond_multiview.png")
    print("\n✨ Disbond region clearly tagged and highlighted in red!")


if __name__ == "__main__":
    main()
