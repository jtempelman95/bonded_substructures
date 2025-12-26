"""Mesh visualization utilities using pyvista and matplotlib."""

from pathlib import Path
from typing import Optional, Union

import gmsh
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from bonded_substructures.mesh.markers import get_tag_name


def plot_mesh_pyvista(
    mesh_file: Optional[Union[str, Path]] = None,
    show_edges: bool = True,
    show_labels: bool = False,
    notebook: bool = False,
) -> Optional[pv.Plotter]:
    """Visualize mesh using PyVista.

    Args:
        mesh_file: Path to .msh file. If None, uses current gmsh model
        show_edges: Show element edges
        show_labels: Show physical group labels
        notebook: Use notebook-friendly plotting

    Returns:
        PyVista plotter object (if notebook=True, returns None after showing)
    """
    if mesh_file is not None:
        # Read from file
        mesh_file = Path(mesh_file)
        if not mesh_file.exists():
            raise FileNotFoundError(f"Mesh file not found: {mesh_file}")

        # Read mesh with pyvista
        mesh = pv.read(str(mesh_file))
    else:
        # Get mesh from gmsh
        if not gmsh.isInitialized():
            raise RuntimeError("Gmsh not initialized and no mesh file provided")

        # Write temporary file and read it
        temp_file = Path("temp_mesh.msh")
        gmsh.write(str(temp_file))
        mesh = pv.read(str(temp_file))
        temp_file.unlink()  # Delete temporary file

    # Create plotter
    plotter = pv.Plotter(notebook=notebook)

    # Plot mesh
    plotter.add_mesh(
        mesh,
        show_edges=show_edges,
        scalars="gmsh:physical" if "gmsh:physical" in mesh.array_names else None,
        cmap="Set3",
        edge_color="black",
        line_width=0.5,
        opacity=0.8,
    )

    plotter.add_axes()
    plotter.show_grid()

    # Set camera position for better 3D view
    plotter.camera_position = 'isometric'
    plotter.camera.zoom(1.2)

    if notebook:
        plotter.show()
        return None
    else:
        return plotter


def plot_mesh_matplotlib(
    mesh_file: Optional[Union[str, Path]] = None,
    show_physical_groups: bool = True,
    figsize: tuple = (10, 6),
    save_path: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """Visualize mesh using Matplotlib (surface view for 3D meshes).

    Args:
        mesh_file: Path to .msh file. If None, uses current gmsh model
        show_physical_groups: Color elements by physical group
        figsize: Figure size (width, height)
        save_path: Path to save figure. If None, displays interactively

    Returns:
        Matplotlib figure object
    """
    # Get mesh data from gmsh
    if mesh_file is not None:
        gmsh.initialize()
        gmsh.open(str(mesh_file))

    if not gmsh.isInitialized():
        raise RuntimeError("Gmsh not initialized and no mesh file provided")

    # Get nodes
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    node_coords = node_coords.reshape(-1, 3)

    # Create node tag to index mapping
    node_map = {tag: i for i, tag in enumerate(node_tags)}

    # Check if we have 3D or 2D mesh
    has_3d = len(gmsh.model.mesh.getElements(dim=3)[0]) > 0
    mesh_dim = 3 if has_3d else 2

    # Get surface elements for visualization
    element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements(dim=2)

    # Create figure with 3D axes if mesh is 3D
    if mesh_dim == 3:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
    else:
        fig, ax = plt.subplots(figsize=figsize)

    # Plot elements
    for elem_type, elem_tags, elem_nodes in zip(element_types, element_tags, element_node_tags):
        # Get nodes per element
        if elem_type == 2:  # Triangle
            nodes_per_elem = 3
        elif elem_type == 3:  # Quadrangle
            nodes_per_elem = 4
        else:
            continue  # Skip other element types

        # Reshape element nodes
        elem_nodes = elem_nodes.reshape(-1, nodes_per_elem)

        # Plot each element
        for elem_node in elem_nodes:
            # Get node coordinates
            coords = np.array([node_coords[node_map[node]] for node in elem_node])

            # Close the polygon
            coords = np.vstack([coords, coords[0]])

            # Plot element edges
            if mesh_dim == 3:
                ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], "k-", linewidth=0.5, alpha=0.3)
            else:
                ax.plot(coords[:, 0], coords[:, 1], "k-", linewidth=0.5, alpha=0.5)

    # Color elements based on their y-coordinate (height) to distinguish materials
    if show_physical_groups and mesh_dim == 3:
        # Get all surface elements
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        # Find the interface height (estimate from node coordinates)
        y_coords = node_coords[:, 1]
        y_min, y_max = y_coords.min(), y_coords.max()
        y_interface = (y_min + y_max) / 2  # Approximate interface location
        interface_tolerance = (y_max - y_min) * 0.15  # 15% tolerance for interface region

        # Define colors for materials
        color_material_1 = '#FF6B6B'  # Red for bottom material
        color_material_2 = '#4ECDC4'  # Cyan for top material
        color_interface = '#FFE66D'    # Yellow for interface region

        # Plot all surface elements colored by their y-position
        for elem_type, elem_tags, elem_nodes in zip(element_types, element_tags, element_node_tags):
            if elem_type != 2:  # Only process triangles
                continue

            nodes_per_elem = 3
            elem_nodes = elem_nodes.reshape(-1, nodes_per_elem)

            for elem_node in elem_nodes:
                coords = np.array([node_coords[node_map[node]] for node in elem_node])

                # Determine color based on centroid y-coordinate
                centroid_y = coords[:, 1].mean()

                if abs(centroid_y - y_interface) < interface_tolerance:
                    color = color_interface
                elif centroid_y < y_interface:
                    color = color_material_1
                else:
                    color = color_material_2

                # Plot triangular face
                tri = Poly3DCollection([coords], alpha=0.7, facecolor=color, edgecolor='k', linewidth=0.3)
                ax.add_collection3d(tri)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor=color_material_1, label='Material 1 (Bottom)'),
            Patch(facecolor=color_material_2, label='Material 2 (Top)'),
            Patch(facecolor=color_interface, label='Interface/Disbond')
        ]
        ax.legend(handles=legend_elements, loc="upper right")

    elif show_physical_groups and mesh_dim == 2:
        # Original 2D coloring logic
        physical_groups = gmsh.model.getPhysicalGroups(dim=2)
        colors = plt.cm.Set3(np.linspace(0, 1, len(physical_groups)))

        for i, (dim, tag) in enumerate(physical_groups):
            try:
                elem_tags = gmsh.model.mesh.getElements(dim=2, tag=tag)[1]
                if len(elem_tags) == 0:
                    continue

                elem_types, _, elem_nodes = gmsh.model.mesh.getElements(dim=2, tag=tag)
            except Exception as e:
                print(f"Warning: Could not get elements for physical group {tag}: {e}")
                continue

            for elem_type, elem_node_tags in zip(elem_types, elem_nodes):
                if elem_type == 2:
                    nodes_per_elem = 3
                elif elem_type == 3:
                    nodes_per_elem = 4
                else:
                    continue

                elem_node_tags = elem_node_tags.reshape(-1, nodes_per_elem)

                for elem_node in elem_node_tags:
                    coords = np.array([node_coords[node_map[node]] for node in elem_node])
                    ax.fill(coords[:, 0], coords[:, 1], color=colors[i], alpha=0.6)

            name = gmsh.model.getPhysicalName(dim, tag)
            ax.plot([], [], color=colors[i], label=name, linewidth=5)

        ax.legend(loc="upper right")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if mesh_dim == 3:
        ax.set_zlabel("z")
        ax.set_box_aspect([1, 1, 1])
    else:
        ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)
    ax.set_title("Mesh Visualization")

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Figure saved to {save_path}")

    if mesh_file is not None:
        gmsh.finalize()

    return fig


def print_mesh_info():
    """Print information about the current mesh."""
    if not gmsh.isInitialized():
        print("Gmsh not initialized")
        return

    # Get mesh statistics
    node_tags, _, _ = gmsh.model.mesh.getNodes()
    num_nodes = len(node_tags)

    element_types, element_tags, _ = gmsh.model.mesh.getElements()
    num_elements = sum(len(tags) for tags in element_tags)

    # Get physical groups
    physical_groups = gmsh.model.getPhysicalGroups()

    print("=" * 50)
    print("Mesh Information")
    print("=" * 50)
    print(f"Number of nodes: {num_nodes}")
    print(f"Number of elements: {num_elements}")
    print(f"\nPhysical Groups:")

    for dim, tag in physical_groups:
        name = gmsh.model.getPhysicalName(dim, tag)
        dim_name = ["point", "curve", "surface", "volume"][dim]
        print(f"  [{dim_name}] {name} (tag={tag})")

    print("=" * 50)
