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
    )

    plotter.add_axes()
    plotter.show_grid()
    plotter.view_xy()

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
    """Visualize 2D mesh using Matplotlib.

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

    # Get elements (2D elements - triangles and quads)
    element_types, element_tags, element_node_tags = gmsh.model.mesh.getElements(dim=2)

    # Create figure
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
            ax.plot(coords[:, 0], coords[:, 1], "k-", linewidth=0.5, alpha=0.5)

    # Get physical groups if requested
    if show_physical_groups:
        physical_groups = gmsh.model.getPhysicalGroups(dim=2)

        # Define colors for physical groups
        colors = plt.cm.Set3(np.linspace(0, 1, len(physical_groups)))

        for i, (dim, tag) in enumerate(physical_groups):
            # Get elements in this physical group
            elem_tags = gmsh.model.mesh.getElements(dim=2, tag=tag)[1]

            if len(elem_tags) == 0:
                continue

            # Get element types and node tags
            elem_types, _, elem_nodes = gmsh.model.mesh.getElements(dim=2, tag=tag)

            for elem_type, elem_node_tags in zip(elem_types, elem_nodes):
                # Get nodes per element
                if elem_type == 2:  # Triangle
                    nodes_per_elem = 3
                elif elem_type == 3:  # Quadrangle
                    nodes_per_elem = 4
                else:
                    continue

                # Reshape element nodes
                elem_node_tags = elem_node_tags.reshape(-1, nodes_per_elem)

                # Plot each element with color
                for elem_node in elem_node_tags:
                    coords = np.array([node_coords[node_map[node]] for node in elem_node])

                    # Fill element
                    ax.fill(coords[:, 0], coords[:, 1], color=colors[i], alpha=0.6)

            # Add label
            name = gmsh.model.getPhysicalName(dim, tag)
            ax.plot([], [], color=colors[i], label=name, linewidth=5)

        ax.legend(loc="upper right")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
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
