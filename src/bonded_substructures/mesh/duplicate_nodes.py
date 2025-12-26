"""Utility to create non-conformal meshes for contact mechanics.

This module provides functions to duplicate nodes at disbond interfaces,
creating separate node sets for Craig-Bampton substructuring.
"""

import numpy as np
import meshio


def duplicate_disbond_interface_nodes(
    input_mesh_file: str,
    output_mesh_file: str,
    disbond_physical_tag: int = 30,
    material1_physical_tag: int = 1,
    material2_physical_tag: int = 2,
) -> dict:
    """Duplicate nodes at disbond interface to create non-conformal mesh.

    Args:
        input_mesh_file: Path to input .msh file (conformal mesh)
        output_mesh_file: Path to output .msh file (non-conformal mesh)
        disbond_physical_tag: Physical tag for disbond region
        material1_physical_tag: Physical tag for material 1
        material2_physical_tag: Physical tag for material 2

    Returns:
        Dictionary with duplication statistics
    """
    # Read the mesh
    mesh = meshio.read(input_mesh_file)

    # meshio stores cells and cell_data as lists of blocks
    # Find tetrahedral blocks and their physical tags
    tetra_blocks = []
    phys_tags_blocks = []

    for i, cell_block in enumerate(mesh.cells):
        if cell_block.type == "tetra":
            tetra_blocks.append(cell_block.data)
            if "gmsh:physical" in mesh.cell_data:
                phys_tags_blocks.append(mesh.cell_data["gmsh:physical"][i])
            else:
                phys_tags_blocks.append(np.zeros(len(cell_block.data), dtype=int))

    if not tetra_blocks:
        raise ValueError("Mesh must contain tetrahedral elements")

    # Concatenate all tetra blocks
    all_cells = np.vstack(tetra_blocks)
    all_phys_tags = np.concatenate(phys_tags_blocks)

    # Identify elements belonging to each region
    disbond_mask = all_phys_tags == disbond_physical_tag
    mat1_mask = all_phys_tags == material1_physical_tag
    mat2_mask = all_phys_tags == material2_physical_tag

    disbond_elements = all_cells[disbond_mask]
    mat1_elements = all_cells[mat1_mask]
    mat2_elements = all_cells[mat2_mask]

    # Find nodes on disbond boundaries
    disbond_nodes = set(disbond_elements.flatten())

    # Find nodes shared between disbond and each material
    mat1_nodes = set(mat1_elements.flatten())
    mat2_nodes = set(mat2_elements.flatten())

    interface1_nodes = disbond_nodes.intersection(mat1_nodes)  # Disbond-Material1 interface
    interface2_nodes = disbond_nodes.intersection(mat2_nodes)  # Disbond-Material2 interface

    print(f"Found {len(interface1_nodes)} nodes at Material1-Disbond interface")
    print(f"Found {len(interface2_nodes)} nodes at Material2-Disbond interface")

    # Create duplicate nodes
    points = mesh.points.copy()
    n_original_points = len(points)

    # Create mapping from original nodes to duplicates
    node_map_mat1 = {}  # For material 1 elements
    node_map_mat2 = {}  # For material 2 elements

    # Duplicate nodes at interface 1 (for material 1)
    for node in interface1_nodes:
        new_node_id = len(points)
        points = np.vstack([points, mesh.points[node]])
        node_map_mat1[node] = new_node_id

    # Duplicate nodes at interface 2 (for material 2)
    for node in interface2_nodes:
        new_node_id = len(points)
        points = np.vstack([points, mesh.points[node]])
        node_map_mat2[node] = new_node_id

    # Update element connectivity
    new_cells = all_cells.copy()

    # Update material 1 elements to use their duplicate nodes
    for i in np.where(mat1_mask)[0]:
        for j in range(4):  # 4 nodes per tetrahedron
            node = new_cells[i, j]
            if node in node_map_mat1:
                new_cells[i, j] = node_map_mat1[node]

    # Update material 2 elements to use their duplicate nodes
    for i in np.where(mat2_mask)[0]:
        for j in range(4):
            node = new_cells[i, j]
            if node in node_map_mat2:
                new_cells[i, j] = node_map_mat2[node]

    # Preserve original triangle surface elements
    triangle_blocks = []
    triangle_phys_tags = []
    for i, cell_block in enumerate(mesh.cells):
        if cell_block.type == "triangle":
            triangle_blocks.append(("triangle", cell_block.data))
            if "gmsh:physical" in mesh.cell_data:
                triangle_phys_tags.append(mesh.cell_data["gmsh:physical"][i])

    # Create new mesh with updated connectivity
    new_cells_list = [("tetra", new_cells)]
    new_cell_data = {"gmsh:physical": [all_phys_tags]}

    # Add triangles back
    if triangle_blocks:
        new_cells_list = triangle_blocks + [("tetra", new_cells)]
        new_cell_data = {"gmsh:physical": triangle_phys_tags + [all_phys_tags]}

    new_mesh = meshio.Mesh(
        points=points,
        cells=new_cells_list,
        cell_data=new_cell_data,
    )

    # Write output mesh
    meshio.write(output_mesh_file, new_mesh)

    stats = {
        "original_nodes": n_original_points,
        "final_nodes": len(points),
        "duplicated_nodes": len(points) - n_original_points,
        "interface1_nodes": len(interface1_nodes),
        "interface2_nodes": len(interface2_nodes),
    }

    print(f"\n✓ Created non-conformal mesh:")
    print(f"  Original nodes: {stats['original_nodes']}")
    print(f"  Final nodes: {stats['final_nodes']}")
    print(f"  Duplicated: {stats['duplicated_nodes']} nodes")
    print(f"  Saved to: {output_mesh_file}")

    return stats
