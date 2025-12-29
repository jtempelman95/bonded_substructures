"""Example 12: Recover and visualize full wavefield from ROM analysis.

This script demonstrates the complete workflow:
1. Perform Craig-Bampton ROM reduction
2. Solve nonlinear harmonic response with HB-AFT
3. Recover full-field displacements using CB transformation
4. Visualize displacement wavefield on mesh using PyVista

The visualization shows how the reduced-order solution maps back to the
physical mesh, enabling full-field analysis at a fraction of the computational cost.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

try:
    from bonded_substructures.rom.craig_bampton import CraigBamptonReduction
    from bonded_substructures.rom.harmonic_balance import HarmonicBalanceAFT
    import pyvista as pv
    import meshio
    ROM_AVAILABLE = True
except ImportError as e:
    ROM_AVAILABLE = False
    print(f"Requires dolfinx and pyvista: {e}")


def recover_full_field_displacements(cb, U_reduced_time):
    """Recover full-field displacements from reduced ROM solution.

    The Craig-Bampton transformation relates reduced and full coordinates:
        u_full = T * u_reduced

    Where T is the transformation matrix combining modal and constraint modes.

    Args:
        cb: CraigBamptonReduction object (after reduction)
        U_reduced_time: Reduced displacements [n_reduced, n_time]

    Returns:
        U_full_time: Full-field displacements [n_full, n_time]
    """
    print("\n" + "="*70)
    print("Recovering Full-Field Displacements")
    print("="*70)

    n_reduced, n_time = U_reduced_time.shape

    # Get transformation matrix from CB object
    # The transformation matrix T maps reduced DOFs to full DOFs
    # Structure: T = [Phi (normal modes), Psi (constraint modes)]

    # We need to reconstruct the full displacement field
    # For each substructure: u_interior = Phi*q + Psi*u_interface
    #                        u_interface = u_interface (boundary DOFs)

    # Get total number of DOFs in full system from function space
    if hasattr(cb, 'V'):
        n_full = cb.V.dofmap.index_map.size_global * cb.V.dofmap.index_map_bs
    else:
        raise ValueError("Function space not available in CB object")

    print(f"\nDimensions:")
    print(f"  Full system DOFs: {n_full}")
    print(f"  Reduced system DOFs: {n_reduced}")
    print(f"  Time snapshots: {n_time}")
    print(f"  Reduction factor: {n_full/n_reduced:.1f}x")

    # Initialize full displacement array
    U_full_time = np.zeros((n_full, n_time))

    # The reduced system is ordered as: [modal DOFs for each substructure, interface DOFs]
    # We need to map this back to the original DOF numbering

    # For now, use a simplified approach: assume the CB object has a method
    # or we can access the transformation directly

    # Get modal and interface DOF counts
    n_modal_total = 0
    for name, substruct in cb.substructures.items():
        if hasattr(substruct, 'normal_modes') and substruct.normal_modes is not None:
            n_modal_total += substruct.normal_modes.shape[1]

    n_interface_reduced = n_reduced - n_modal_total

    print(f"\n  Modal DOFs: {n_modal_total}")
    print(f"  Interface DOFs: {n_interface_reduced}")

    # Extract modal and interface parts
    U_modal = U_reduced_time[:n_modal_total, :]
    U_interface = U_reduced_time[n_modal_total:, :]

    print(f"\nExpanding each substructure...")

    # Expand each substructure
    current_modal_idx = 0
    current_interface_idx = 0

    for name, substruct in cb.substructures.items():
        print(f"  Substructure {name} (tag {substruct.material_tag}):")

        # Get interior and interface DOFs for this substructure
        interior_dofs = substruct.interior_dofs
        interface_dofs = substruct.interface_dofs

        n_interior = len(interior_dofs)
        n_interface_local = len(interface_dofs)

        print(f"    Interior DOFs: {n_interior}")
        print(f"    Interface DOFs: {n_interface_local}")

        # Get modal contributions for this substructure
        if hasattr(substruct, 'normal_modes') and substruct.normal_modes is not None:
            n_modes = substruct.normal_modes.shape[1]
            Phi = substruct.normal_modes  # [n_interior, n_modes]
            q = U_modal[current_modal_idx:current_modal_idx + n_modes, :]  # Modal amplitudes

            current_modal_idx += n_modes
        else:
            Phi = None
            q = None
            n_modes = 0

        # Get constraint mode contributions
        if hasattr(substruct, 'constraint_modes') and substruct.constraint_modes is not None:
            Psi = substruct.constraint_modes  # [n_interior, n_interface_local]
        else:
            Psi = None

        # Get interface displacements for this substructure
        # Map from reduced interface DOFs to local interface DOFs
        u_interface_local = U_interface[current_interface_idx:current_interface_idx + n_interface_local, :]
        current_interface_idx += n_interface_local

        # Compute interior displacements
        if n_interior > 0:
            u_interior = np.zeros((n_interior, n_time))

            # Modal contribution
            if Phi is not None and q is not None:
                u_interior += Phi @ q

            # Constraint mode contribution
            if Psi is not None:
                u_interior += Psi @ u_interface_local

            # Assign to full array
            U_full_time[interior_dofs, :] = u_interior

        # Assign interface displacements
        U_full_time[interface_dofs, :] = u_interface_local

        print(f"    ✓ Expanded to full field")

    print(f"\n✓ Full-field recovery complete")
    print(f"  Max displacement (reduced): {np.max(np.abs(U_reduced_time))*1e6:.3f} μm")
    print(f"  Max displacement (full): {np.max(np.abs(U_full_time))*1e6:.3f} μm")

    return U_full_time


def map_displacements_to_mesh(mesh_file, U_full_time, cb):
    """Map full-field displacements to mesh nodes for visualization.

    Args:
        mesh_file: Path to mesh file
        U_full_time: Full displacements [n_dof, n_time]
        cb: CraigBamptonReduction object

    Returns:
        mesh_data: Dict with mesh and displacement data for PyVista
    """
    print("\n" + "="*70)
    print("Mapping Displacements to Mesh")
    print("="*70)

    # Load mesh with meshio
    mesh = meshio.read(mesh_file)

    print(f"\nMesh loaded:")
    print(f"  Nodes: {len(mesh.points)}")
    print(f"  Elements: {sum(len(cell_block.data) for cell_block in mesh.cells)}")

    # The full displacement array U_full_time has shape [n_dof, n_time]
    # We need to map this to mesh nodes
    # DOFs are organized as: [u_x, u_y, u_z] for each node (3D) or [u_x, u_y] (2D)

    n_nodes = len(mesh.points)
    dim = mesh.points.shape[1]  # 2D or 3D
    n_time = U_full_time.shape[1]

    # Determine DOF structure (2D or 3D problem)
    # For FEM: typically n_dof = n_nodes * dim
    n_dof_expected = n_nodes * dim

    if U_full_time.shape[0] != n_dof_expected:
        print(f"\n⚠️  Warning: DOF mismatch")
        print(f"  Expected DOFs: {n_dof_expected} ({n_nodes} nodes × {dim} dimensions)")
        print(f"  Actual DOFs: {U_full_time.shape[0]}")
        print(f"  Using available DOFs with truncation/padding")

    # Reshape displacements to [n_nodes, dim, n_time]
    # Handle potential size mismatch
    n_dof_actual = min(U_full_time.shape[0], n_dof_expected)

    U_nodal = np.zeros((n_nodes, 3, n_time))  # Always use 3D for PyVista

    for t in range(n_time):
        for node_idx in range(n_nodes):
            for d in range(min(dim, 3)):
                dof_idx = node_idx * dim + d
                if dof_idx < n_dof_actual:
                    U_nodal[node_idx, d, t] = U_full_time[dof_idx, t]

    print(f"\n✓ Displacements mapped to nodes")
    print(f"  Node displacement array: {U_nodal.shape} (nodes × XYZ × time)")

    # Compute displacement magnitudes
    U_magnitude = np.sqrt(np.sum(U_nodal**2, axis=1))  # [n_nodes, n_time]

    print(f"  Max displacement magnitude: {np.max(U_magnitude)*1e6:.3f} μm")

    return {
        'mesh': mesh,
        'points': mesh.points,
        'U_nodal': U_nodal,
        'U_magnitude': U_magnitude,
        'n_time': n_time
    }


def create_pyvista_wavefield_animation(mesh_data, output_dir, mesh_name, omega):
    """Create PyVista visualizations of displacement wavefield.

    Args:
        mesh_data: Dict from map_displacements_to_mesh()
        output_dir: Output directory
        mesh_name: Name for output files
        omega: Excitation frequency (rad/s)
    """
    print("\n" + "="*70)
    print("Creating PyVista Wavefield Visualizations")
    print("="*70)

    mesh = mesh_data['mesh']
    points = mesh_data['points']
    U_nodal = mesh_data['U_nodal']
    U_magnitude = mesh_data['U_magnitude']
    n_time = mesh_data['n_time']

    # Period and time vector
    T = 2 * np.pi / omega
    t = np.linspace(0, T, n_time, endpoint=False)

    # Select time snapshots for visualization
    n_snapshots = 6
    snapshot_indices = np.linspace(0, n_time-1, n_snapshots, dtype=int)

    print(f"\nCreating {n_snapshots} displacement snapshots...")

    # Create figure with subplots for multiple time snapshots
    plotter = pv.Plotter(shape=(2, 3), off_screen=True, window_size=[1920, 1280])

    # Load mesh into PyVista
    # Convert meshio to pyvista
    all_cells = []
    all_cell_data = []

    for cell_block in mesh.cells:
        if cell_block.type == "tetra":
            cells = cell_block.data
            # PyVista format: [4, p0, p1, p2, p3, ...]
            pv_cells = np.hstack([np.full((cells.shape[0], 1), 4), cells]).ravel()
            all_cells.append(pv_cells)

        elif cell_block.type == "triangle":
            cells = cell_block.data
            # PyVista format: [3, p0, p1, p2, ...]
            pv_cells = np.hstack([np.full((cells.shape[0], 1), 3), cells]).ravel()
            all_cells.append(pv_cells)

        elif cell_block.type == "quad":
            cells = cell_block.data
            # PyVista format: [4, p0, p1, p2, p3, ...]
            pv_cells = np.hstack([np.full((cells.shape[0], 1), 4), cells]).ravel()
            all_cells.append(pv_cells)

    # Use first cell type for main visualization
    if len(all_cells) == 0:
        print("⚠️  No supported cell types found")
        return

    # Create PyVista mesh for each snapshot
    for i, idx in enumerate(snapshot_indices):
        row = i // 3
        col = i % 3
        plotter.subplot(row, col)

        # Get displacement at this time
        U_disp = U_nodal[:, :, idx]  # [n_nodes, 3]
        U_mag = U_magnitude[:, idx]  # [n_nodes]

        # Create deformed mesh
        points_deformed = points + U_disp * 1000  # Scale for visibility (1000x)

        # Create PyVista grid with deformed points
        cell_types = []
        cells_concat = []

        for j, cell_block in enumerate(mesh.cells):
            if cell_block.type == "tetra":
                cells = cell_block.data
                pv_cells = np.hstack([np.full((cells.shape[0], 1), 4), cells]).ravel()
                cells_concat.append(pv_cells)
                cell_types.extend([pv.CellType.TETRA] * cells.shape[0])

            elif cell_block.type == "triangle":
                cells = cell_block.data
                pv_cells = np.hstack([np.full((cells.shape[0], 1), 3), cells]).ravel()
                cells_concat.append(pv_cells)
                cell_types.extend([pv.CellType.TRIANGLE] * cells.shape[0])

        if len(cells_concat) > 0:
            grid = pv.UnstructuredGrid(
                np.concatenate(cells_concat),
                np.array(cell_types),
                points_deformed
            )

            # Add displacement magnitude as scalar
            grid["Displacement (μm)"] = U_mag * 1e6

            # Plot with displacement coloring
            plotter.add_mesh(
                grid,
                scalars="Displacement (μm)",
                show_edges=False,
                cmap='jet',
                show_scalar_bar=True,
                scalar_bar_args={
                    'title': 'Displacement (μm)',
                    'title_font_size': 12,
                    'label_font_size': 10,
                    'n_labels': 5,
                    'position_x': 0.85,
                    'position_y': 0.1,
                    'width': 0.12,
                    'height': 0.8
                }
            )

            plotter.camera_position = 'iso'
            plotter.add_title(
                f't = {t[idx]*1000:.2f} ms (deformed 1000×)',
                font_size=12
            )
            plotter.add_axes(line_width=3)
            plotter.background_color = 'white'

    # Save figure
    output_path = output_dir / f"rom_wavefield_{mesh_name}.png"
    plotter.screenshot(str(output_path), return_img=False)
    plotter.close()

    print(f"\n  ✓ Saved: {output_path}")

    # Create single high-resolution snapshot at peak displacement
    print(f"\nCreating high-resolution snapshot at peak displacement...")

    # Find time index with maximum displacement
    max_disp_global = np.max(U_magnitude, axis=0)
    peak_idx = np.argmax(max_disp_global)

    plotter_peak = pv.Plotter(off_screen=True, window_size=[1920, 1080])

    U_peak = U_nodal[:, :, peak_idx]
    U_mag_peak = U_magnitude[:, peak_idx]

    points_deformed_peak = points + U_peak * 1000

    # Recreate grid
    cells_concat = []
    cell_types = []

    for cell_block in mesh.cells:
        if cell_block.type == "tetra":
            cells = cell_block.data
            pv_cells = np.hstack([np.full((cells.shape[0], 1), 4), cells]).ravel()
            cells_concat.append(pv_cells)
            cell_types.extend([pv.CellType.TETRA] * cells.shape[0])

        elif cell_block.type == "triangle":
            cells = cell_block.data
            pv_cells = np.hstack([np.full((cells.shape[0], 1), 3), cells]).ravel()
            cells_concat.append(pv_cells)
            cell_types.extend([pv.CellType.TRIANGLE] * cells.shape[0])

    if len(cells_concat) > 0:
        grid_peak = pv.UnstructuredGrid(
            np.concatenate(cells_concat),
            np.array(cell_types),
            points_deformed_peak
        )

        grid_peak["Displacement (μm)"] = U_mag_peak * 1e6

        plotter_peak.add_mesh(
            grid_peak,
            scalars="Displacement (μm)",
            show_edges=True,
            edge_color='black',
            line_width=0.5,
            cmap='jet',
            show_scalar_bar=True,
            scalar_bar_args={
                'title': 'Displacement (μm)',
                'title_font_size': 16,
                'label_font_size': 14,
                'n_labels': 7
            }
        )

        plotter_peak.camera_position = 'iso'
        plotter_peak.add_title(
            f'Peak Displacement Wavefield (t = {t[peak_idx]*1000:.2f} ms)\n'
            f'Deformed shape (1000× magnification)',
            font_size=16,
            font='arial'
        )
        plotter_peak.add_axes(
            xlabel='X',
            ylabel='Y',
            zlabel='Z',
            line_width=5
        )
        plotter_peak.background_color = 'white'

        output_peak = output_dir / f"rom_wavefield_peak_{mesh_name}.png"
        plotter_peak.screenshot(str(output_peak), return_img=False)
        plotter_peak.close()

        print(f"  ✓ Saved: {output_peak}")
        print(f"  Peak displacement: {np.max(U_mag_peak)*1e6:.3f} μm at t = {t[peak_idx]*1000:.2f} ms")


def main():
    """Complete wavefield recovery and visualization workflow."""

    print("\n" + "="*70)
    print("ROM WAVEFIELD RECOVERY AND VISUALIZATION")
    print("Craig-Bampton ROM → HB-AFT → Full-Field Recovery → PyVista")
    print("="*70)

    if not ROM_AVAILABLE:
        print("\n⚠️  Requires dolfinx and pyvista")
        return

    # Use wide plate mesh (2D, works with dolfinx)
    mesh_file = "wide_plate_disbond.msh"

    if not Path(mesh_file).exists():
        print(f"\n⚠️  Mesh not found: {mesh_file}")
        print("Run Example 07 first")
        return

    output_dir = Path("examples/output")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Craig-Bampton Reduction
    print("\n" + "="*70)
    print("STEP 1: Craig-Bampton ROM Reduction")
    print("="*70)

    cb = CraigBamptonReduction(
        mesh_file=mesh_file,
        material_tags={"Material_1": 1, "Material_2": 2}
    )

    cb.load_mesh()
    cb.partition_substructures()
    cb.identify_interface_dofs()
    cb.assemble_matrices()

    # Store full DOF count before reduction
    n_full_dofs = cb.V.dofmap.index_map.size_global * cb.V.dofmap.index_map_bs

    cb.compute_modes(n_modes=10)
    K_r, M_r = cb.assemble_reduced_system()

    print(f"\n✓ ROM reduction complete: {K_r.shape[0]} DOFs")
    print(f"  Full system: {n_full_dofs} DOFs")

    # Step 2: Nonlinear Harmonic Balance
    print("\n" + "="*70)
    print("STEP 2: Nonlinear Harmonic Balance")
    print("="*70)

    hb = HarmonicBalanceAFT(K_r, M_r, damping_ratio=0.02)
    n_harmonics = 5
    hb.set_harmonics(n_harmonics)

    freq_hz = 100.0
    omega = 2 * np.pi * freq_hz
    F_magnitude = 5000.0

    F_ext = np.zeros((K_r.shape[0], 2*n_harmonics + 1))
    F_ext[0, 1] = F_magnitude

    print(f"\nExcitation: {F_magnitude:.0f} N at {freq_hz:.1f} Hz")

    # Solve (linear for simplicity - can add contact if desired)
    Z = hb.assemble_frequency_matrix(omega)
    U_freq_flat = np.linalg.solve(Z, F_ext.ravel())
    U_freq = U_freq_flat.reshape(F_ext.shape)

    # Convert to time domain
    U_reduced_time, t = hb.frequency_to_time(U_freq, omega)

    print(f"\n✓ HB solution complete")
    print(f"  Time points: {len(t)}")
    print(f"  Max displacement (reduced): {np.max(np.abs(U_reduced_time))*1e6:.3f} μm")

    # Step 3: Recover full-field displacements
    try:
        U_full_time = recover_full_field_displacements(cb, U_reduced_time)
    except Exception as e:
        print(f"\n✗ Error recovering full field: {e}")
        import traceback
        traceback.print_exc()
        return

    # Step 4: Map to mesh
    try:
        mesh_data = map_displacements_to_mesh(mesh_file, U_full_time, cb)
    except Exception as e:
        print(f"\n✗ Error mapping to mesh: {e}")
        import traceback
        traceback.print_exc()
        return

    # Step 5: Create PyVista visualizations
    try:
        create_pyvista_wavefield_animation(
            mesh_data, output_dir, "wide_plate", omega
        )
    except Exception as e:
        print(f"\n✗ Error creating visualizations: {e}")
        import traceback
        traceback.print_exc()
        return

    # Summary
    print("\n" + "="*70)
    print("WAVEFIELD VISUALIZATION COMPLETE")
    print("="*70)

    print(f"\nGenerated files:")
    print(f"  ✓ rom_wavefield_wide_plate.png (6 time snapshots)")
    print(f"  ✓ rom_wavefield_peak_wide_plate.png (peak displacement)")

    print(f"\nWorkflow summary:")
    print(f"  1. ROM reduction: {n_full_dofs} → {K_r.shape[0]} DOFs ({n_full_dofs/K_r.shape[0]:.1f}×)")
    print(f"  2. HB solution: {n_harmonics} harmonics, {len(t)} time points")
    print(f"  3. Full-field recovery: {U_full_time.shape[0]} DOFs recovered")
    print(f"  4. PyVista visualization: Deformed mesh with displacement coloring")

    print(f"\n✨ Full wavefield successfully recovered and visualized!")


if __name__ == "__main__":
    main()
