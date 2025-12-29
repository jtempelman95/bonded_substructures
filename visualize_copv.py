"""Create comprehensive visualizations of ideal and disbonded COPV meshes."""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

try:
    import pyvista as pv
    import meshio
    PYVISTA_AVAILABLE = True
except ImportError:
    PYVISTA_AVAILABLE = False
    print("PyVista not available - will create matplotlib-only visualizations")


def load_mesh_from_gmsh(mesh_file):
    """Load mesh from gmsh .msh file using meshio and convert to pyvista."""
    mesh = meshio.read(mesh_file)

    # Extract points and cells
    points = mesh.points

    # Get tetrahedra
    cells = []
    for cell_block in mesh.cells:
        if cell_block.type == "tetra":
            cells = cell_block.data
            break

    # Get cell data (physical tags)
    cell_data = {}
    if "gmsh:physical" in mesh.cell_data:
        for i, cell_block in enumerate(mesh.cells):
            if cell_block.type == "tetra":
                cell_data["Material"] = mesh.cell_data["gmsh:physical"][i]
                break

    return points, cells, cell_data


def create_pyvista_visualization(mesh_file, title, output_file):
    """Create PyVista visualization of COPV mesh."""
    if not PYVISTA_AVAILABLE:
        print(f"  Skipping PyVista visualization (not available)")
        return

    # Read mesh with meshio
    mesh = meshio.read(mesh_file)

    # Convert to pyvista
    points = mesh.points
    cells_dict = {}
    for cell_block in mesh.cells:
        cells_dict[cell_block.type] = cell_block.data

    # Create pyvista mesh
    if "tetra" in cells_dict:
        cells = cells_dict["tetra"]
        # PyVista format: [n_points, p0, p1, p2, p3, n_points, ...]
        pv_cells = np.hstack([np.full((cells.shape[0], 1), 4), cells]).ravel()
        pv_mesh = pv.UnstructuredGrid(pv_cells, np.full(cells.shape[0], 10), points)

        # Add material data if available
        if "gmsh:physical" in mesh.cell_data:
            for i, cell_block in enumerate(mesh.cells):
                if cell_block.type == "tetra":
                    pv_mesh["Material"] = mesh.cell_data["gmsh:physical"][i]
                    break
    else:
        print(f"  Warning: No tetrahedral cells found in {mesh_file}")
        return

    # Create plotter
    pl = pv.Plotter(off_screen=True, window_size=[1920, 1080])

    # Add mesh with material coloring
    if "Material" in pv_mesh.array_names:
        pl.add_mesh(pv_mesh, scalars="Material", show_edges=True,
                   edge_color="black", line_width=0.5,
                   cmap="viridis", show_scalar_bar=True,
                   scalar_bar_args={'title': 'Material Tag'})
    else:
        pl.add_mesh(pv_mesh, color="lightblue", show_edges=True, edge_color="black")

    # Set camera view
    pl.camera_position = 'iso'
    pl.add_title(title, font_size=16)
    pl.add_axes()

    # Save screenshot
    pl.screenshot(output_file, return_img=False)
    pl.close()

    print(f"  ✓ Saved: {output_file}")


def create_matplotlib_cross_section(mesh_file, title, output_file):
    """Create matplotlib cross-section visualization."""
    # Load mesh
    points, cells, cell_data = load_mesh_from_gmsh(mesh_file)

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle(title, fontsize=16, fontweight='bold')

    # Extract coordinates
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]

    # Plot 1: XY plane (top view)
    ax = axes[0, 0]
    scatter = ax.scatter(x, y, c=z, cmap='viridis', s=1, alpha=0.5)
    ax.set_xlabel('X (m)', fontsize=11)
    ax.set_ylabel('Y (m)', fontsize=11)
    ax.set_title('Top View (XY Plane)', fontsize=12, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='Z (m)')

    # Plot 2: XZ plane (side view)
    ax = axes[0, 1]
    scatter = ax.scatter(x, z, c=y, cmap='plasma', s=1, alpha=0.5)
    ax.set_xlabel('X (m)', fontsize=11)
    ax.set_ylabel('Z (m)', fontsize=11)
    ax.set_title('Side View (XZ Plane)', fontsize=12, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='Y (m)')

    # Plot 3: YZ plane (front view)
    ax = axes[1, 0]
    scatter = ax.scatter(y, z, c=x, cmap='coolwarm', s=1, alpha=0.5)
    ax.set_xlabel('Y (m)', fontsize=11)
    ax.set_ylabel('Z (m)', fontsize=11)
    ax.set_title('Front View (YZ Plane)', fontsize=12, fontweight='bold')
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='X (m)')

    # Plot 4: Radial distribution
    ax = axes[1, 1]
    r = np.sqrt(x**2 + y**2)
    scatter = ax.scatter(r, z, c=z, cmap='viridis', s=1, alpha=0.5)
    ax.set_xlabel('Radial Distance (m)', fontsize=11)
    ax.set_ylabel('Z (m)', fontsize=11)
    ax.set_title('Radial Distribution', fontsize=12, fontweight='bold')
    ax.grid(True, alpha=0.3)
    plt.colorbar(scatter, ax=ax, label='Z (m)')

    # Add mesh statistics
    stats_text = f"Nodes: {len(points):,}\nElements: {len(cells):,}"
    fig.text(0.02, 0.02, stats_text, fontsize=10, family='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  ✓ Saved: {output_file}")


def main():
    """Create all COPV visualizations."""

    print("=" * 70)
    print("COPV Mesh Visualization Suite")
    print("=" * 70)

    # Ensure output directory exists
    output_dir = Path("examples/output")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Define meshes to visualize
    meshes = [
        {
            "file": "bonded_cylinder_basic.msh",
            "title": "Ideal COPV (No Disbond)",
            "short": "ideal"
        },
        {
            "file": "bonded_cylinder_disbond.msh",
            "title": "COPV with Circumferential Disbond",
            "short": "disbond"
        }
    ]

    for mesh_info in meshes:
        mesh_file = mesh_info["file"]
        title = mesh_info["title"]
        short_name = mesh_info["short"]

        if not Path(mesh_file).exists():
            print(f"\n⚠️  Mesh file not found: {mesh_file}")
            print(f"   Run examples/04_bonded_cylinder_basic.py or 05_bonded_cylinder_disbond.py first")
            continue

        print(f"\n{title}")
        print("-" * 70)

        # PyVista 3D visualization
        if PYVISTA_AVAILABLE:
            pv_output = output_dir / f"copv_{short_name}_3d.png"
            print(f"Creating 3D visualization...")
            try:
                create_pyvista_visualization(mesh_file, title, str(pv_output))
            except Exception as e:
                print(f"  ✗ PyVista visualization failed: {e}")

        # Matplotlib cross-sections
        mpl_output = output_dir / f"copv_{short_name}_sections.png"
        print(f"Creating cross-section plots...")
        try:
            create_matplotlib_cross_section(mesh_file, title, str(mpl_output))
        except Exception as e:
            print(f"  ✗ Matplotlib visualization failed: {e}")
            import traceback
            traceback.print_exc()

    print("\n" + "=" * 70)
    print("Visualization Complete!")
    print("=" * 70)
    print(f"\nGenerated files in: {output_dir}/")
    print("  - copv_ideal_3d.png (3D mesh view)")
    print("  - copv_ideal_sections.png (cross-sections)")
    print("  - copv_disbond_3d.png (3D mesh view)")
    print("  - copv_disbond_sections.png (cross-sections)")


if __name__ == "__main__":
    main()
