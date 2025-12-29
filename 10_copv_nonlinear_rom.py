"""Example 10: Nonlinear ROM analysis of COPV with contact dynamics.

This example demonstrates the complete ROM workflow for COPV meshes:
1. Load ideal and disbonded COPV meshes
2. Perform Craig-Bampton reduction
3. Run nonlinear harmonic balance with contact using HB-AFT
4. Visualize displacement fields at time snapshots

The analysis captures the nonlinear dynamics of disbond contact in the
circumferential interface region of the composite overwrapped pressure vessel.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

try:
    from bonded_substructures.rom.craig_bampton import CraigBamptonReduction
    from bonded_substructures.rom.harmonic_balance import HarmonicBalanceAFT
    ROM_AVAILABLE = True
except ImportError as e:
    ROM_AVAILABLE = False
    print(f"ROM requires dolfinx: {e}")


def perform_rom_reduction(mesh_file, n_modes=10):
    """Perform Craig-Bampton reduction on COPV mesh.

    Args:
        mesh_file: Path to .msh file
        n_modes: Number of normal modes per substructure

    Returns:
        tuple: (K_r, M_r, cb) - Reduced matrices and CB object
    """
    print(f"\n{'='*70}")
    print(f"Craig-Bampton Reduction: {mesh_file}")
    print('='*70)

    cb = CraigBamptonReduction(
        mesh_file=mesh_file,
        material_tags={"Material_1": 1, "Material_2": 2}
    )

    # Load and partition mesh
    print("\nStep 1: Loading mesh...")
    cb.load_mesh()

    print("Step 2: Partitioning substructures...")
    cb.partition_substructures()

    print("Step 3: Identifying interface DOFs...")
    cb.identify_interface_dofs()

    print("Step 4: Assembling global matrices...")
    cb.assemble_matrices()

    # Compute modes
    print(f"\nStep 5: Computing {n_modes} Craig-Bampton modes per substructure...")
    cb.compute_modes(n_modes=n_modes)

    # Assemble reduced system
    print("\nStep 6: Assembling reduced system...")
    K_r, M_r = cb.assemble_reduced_system()

    # Print reduction statistics
    n_full = cb.K.shape[0] if hasattr(cb, 'K') else 'unknown'
    n_reduced = K_r.shape[0]
    if isinstance(n_full, int):
        reduction_factor = n_full / n_reduced
        print(f"\n{'='*70}")
        print(f"Reduction Complete:")
        print(f"  Original DOFs: {n_full:,}")
        print(f"  Reduced DOFs: {n_reduced}")
        print(f"  Reduction factor: {reduction_factor:.1f}x")
        print('='*70)
    else:
        print(f"\n  Reduced system: {n_reduced} DOFs")

    return K_r, M_r, cb


def run_nonlinear_hb_analysis(K_r, M_r, mesh_name, has_disbond=False):
    """Run nonlinear harmonic balance analysis with contact.

    Args:
        K_r: Reduced stiffness matrix
        M_r: Reduced mass matrix
        mesh_name: Name for output files
        has_disbond: Whether this configuration has disbond

    Returns:
        tuple: (U_freq, u_time, f_nl_time, hb, contact_params, omega) - Results
    """
    print(f"\n{'='*70}")
    print(f"Nonlinear Harmonic Balance Analysis: {mesh_name}")
    print('='*70)

    # Initialize HB solver
    damping_ratio = 0.02  # 2% damping
    hb = HarmonicBalanceAFT(K_r, M_r, damping_ratio=damping_ratio)

    # Setup harmonics
    n_harmonics = 7  # Use 7 harmonics for better accuracy
    hb.set_harmonics(n_harmonics)

    print(f"\nHarmonic Balance Setup:")
    print(f"  Reduced DOFs: {K_r.shape[0]}")
    print(f"  Damping ratio: {damping_ratio*100:.1f}%")
    print(f"  Number of harmonics: {n_harmonics}")
    print(f"  Time points per period: {hb.n_time_points}")
    print(f"  Total unknowns: {K_r.shape[0] * (2*n_harmonics + 1):,}")

    # Excitation parameters - realistic for COPV vibration
    freq_hz = 50.0  # Hz (typical structural frequency for COPV)
    omega = 2 * np.pi * freq_hz
    F_magnitude = 5000.0  # N (moderate harmonic load)

    print(f"\nExcitation:")
    print(f"  Frequency: {freq_hz:.1f} Hz ({omega:.1f} rad/s)")
    print(f"  Force magnitude: {F_magnitude:,.0f} N")
    print(f"  Period: {1000/freq_hz:.2f} ms")

    # Setup external force (harmonic at fundamental frequency)
    F_ext = np.zeros((K_r.shape[0], 2*n_harmonics + 1))
    # Apply force at first few DOFs (boundary excitation)
    for i in range(min(3, K_r.shape[0])):
        F_ext[i, 1] = F_magnitude / 3  # Cosine component, distributed

    # Contact parameters
    if has_disbond:
        # For disbonded case, assume interface DOFs may have contact
        # Estimate: last 20% of DOFs are interface-related
        n_interface_estimate = int(K_r.shape[0] * 0.2)
        start_idx = K_r.shape[0] - n_interface_estimate
        # Select subset for contact (avoid too many for stability)
        contact_dofs = list(range(start_idx, K_r.shape[0], max(1, n_interface_estimate // 10)))[:10]

        contact_params = {
            'contact_dofs': contact_dofs,
            'contact_stiffness': 1e7,  # N/m (moderate stiffness for ROM)
            'gap_initial': 5e-6,  # 5 μm initial gap (realistic for disbond)
            'contact_type': 'penalty'
        }

        print(f"\nContact Parameters:")
        print(f"  Contact DOFs: {len(contact_dofs)} (indices: {contact_dofs[0]}-{contact_dofs[-1]})")
        print(f"  Contact stiffness: {contact_params['contact_stiffness']:.1e} N/m")
        print(f"  Initial gap: {contact_params['gap_initial']*1e6:.1f} μm")

    else:
        # No contact for ideal case (but still run to show baseline)
        contact_params = {
            'contact_dofs': [],
            'contact_stiffness': 0.0,
            'gap_initial': 0.0,
            'contact_type': 'penalty'
        }
        print(f"\n  No contact (ideal COPV)")

    # Solve nonlinear harmonic balance
    print(f"\n{'='*70}")
    print("Solving Nonlinear Harmonic Balance (AFT Method)...")
    print('='*70)

    if has_disbond:
        # Nonlinear solution with contact
        U_freq, u_time, history = hb.solve_harmonic_balance(
            omega, F_ext, contact_params,
            max_iter=150, tol=1e-6, relaxation=0.4
        )

        # Compute contact forces
        f_nl_time = hb.evaluate_contact_forces(u_time, contact_params)

        print(f"\nConvergence:")
        print(f"  Iterations: {len(history['iteration'])}")
        print(f"  Final residual: {history['residual'][-1]:.3e}")

    else:
        # Linear solution (no contact)
        # Build frequency domain system matrix
        Z = hb.assemble_frequency_matrix(omega)

        # Solve linear system
        F_ext_flat = F_ext.ravel()
        U_freq_flat = np.linalg.solve(Z, F_ext_flat)
        U_freq = U_freq_flat.reshape(F_ext.shape)

        # Convert to time domain
        u_time = hb.frequency_to_time(U_freq, omega)

        # No contact forces
        f_nl_time = np.zeros_like(u_time)

        print(f"\n  Linear solution (no iterations needed)")

    # Print response statistics
    print(f"\nResponse Statistics:")
    print(f"  Max displacement: {np.max(np.abs(u_time))*1e6:.3f} μm")
    print(f"  RMS displacement: {np.sqrt(np.mean(u_time**2))*1e6:.3f} μm")

    if has_disbond:
        print(f"  Max contact force: {np.max(f_nl_time):.3e} N")
        contact_active = np.sum(f_nl_time > 1e-6)
        contact_pct = contact_active / f_nl_time.size * 100
        print(f"  Contact active: {contact_pct:.2f}% of time-DOF pairs")

    return U_freq, u_time, f_nl_time, hb, contact_params, omega


def visualize_displacement_snapshots(u_time, hb, omega, mesh_name, output_dir, contact_params=None, f_nl_time=None):
    """Create comprehensive visualization of displacement fields and contact.

    Args:
        u_time: Time-domain displacements [n_dof, n_time]
        hb: HarmonicBalanceAFT object
        omega: Excitation frequency (rad/s)
        mesh_name: Name for plot title
        output_dir: Output directory for plots
        contact_params: Contact parameters dict (optional)
        f_nl_time: Contact forces [n_dof, n_time] (optional)
    """
    print(f"\n{'='*70}")
    print("Creating Displacement Field Visualizations")
    print('='*70)

    # Time vector
    T = 2 * np.pi / omega
    t = np.linspace(0, T, hb.n_time_points, endpoint=False)

    # Select time snapshots: 0, T/4, T/2, 3T/4
    snapshot_indices = [0, hb.n_time_points//4, hb.n_time_points//2, 3*hb.n_time_points//4]
    snapshot_labels = ['t = 0', 't = T/4', 't = T/2', 't = 3T/4']

    # Create figure with multiple subplots
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(4, 3, hspace=0.35, wspace=0.3)

    # Plot 1-4: Displacement field snapshots
    for i, (idx, label) in enumerate(zip(snapshot_indices, snapshot_labels)):
        ax = fig.add_subplot(gs[i//2, i%2])

        # Get displacement at this time
        u_snapshot = u_time[:, idx]
        dof_indices = np.arange(len(u_snapshot))

        # Plot displacement field
        ax.plot(dof_indices, u_snapshot * 1e6, 'b-', linewidth=1.5, alpha=0.8)
        ax.fill_between(dof_indices, 0, u_snapshot * 1e6, alpha=0.3)

        # Highlight contact DOFs if present
        if contact_params and len(contact_params['contact_dofs']) > 0:
            contact_dofs = contact_params['contact_dofs']
            ax.scatter(contact_dofs, u_snapshot[contact_dofs] * 1e6,
                      color='red', s=50, zorder=5, label='Contact DOFs')
            ax.axhline(contact_params['gap_initial'] * 1e6, color='k',
                      linestyle='--', linewidth=1, alpha=0.5, label='Gap threshold')
            ax.legend(fontsize=8, loc='upper right')

        ax.set_xlabel('DOF Index', fontsize=10)
        ax.set_ylabel('Displacement (μm)', fontsize=10)
        ax.set_title(f'Displacement Field: {label}', fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)

    # Plot 5: Time history at selected DOFs
    ax5 = fig.add_subplot(gs[0, 2])
    dof_samples = [0, len(u_time)//4, len(u_time)//2, 3*len(u_time)//4]
    for dof in dof_samples[:4]:
        if dof < len(u_time):
            ax5.plot(t * 1000, u_time[dof, :] * 1e6, linewidth=1.5, alpha=0.8,
                    label=f'DOF {dof}')
    ax5.set_xlabel('Time (ms)', fontsize=10)
    ax5.set_ylabel('Displacement (μm)', fontsize=10)
    ax5.set_title('Time History (Selected DOFs)', fontsize=11, fontweight='bold')
    ax5.legend(fontsize=8)
    ax5.grid(True, alpha=0.3)

    # Plot 6: Contact forces (if available)
    ax6 = fig.add_subplot(gs[1, 2])
    if f_nl_time is not None and contact_params and len(contact_params['contact_dofs']) > 0:
        total_contact = np.sum(f_nl_time[contact_params['contact_dofs'], :], axis=0)
        ax6.plot(t * 1000, total_contact, 'r-', linewidth=2, label='Total Contact Force')
        ax6.fill_between(t * 1000, 0, total_contact, alpha=0.3, color='red')

        # Shade contact active regions
        contact_active = total_contact > 1e-3
        for i in range(len(t)-1):
            if contact_active[i]:
                ax6.axvspan(t[i]*1000, t[i+1]*1000, alpha=0.15, color='red')

        ax6.set_xlabel('Time (ms)', fontsize=10)
        ax6.set_ylabel('Contact Force (N)', fontsize=10)
        ax6.set_title('Total Contact Force', fontsize=11, fontweight='bold')
        ax6.grid(True, alpha=0.3)
    else:
        ax6.text(0.5, 0.5, 'No Contact\n(Ideal COPV)',
                horizontalalignment='center', verticalalignment='center',
                transform=ax6.transAxes, fontsize=14, fontweight='bold')
        ax6.set_title('Contact Force', fontsize=11, fontweight='bold')

    # Plot 7: Maximum displacement over time
    ax7 = fig.add_subplot(gs[2, 2])
    max_disp = np.max(np.abs(u_time), axis=0)
    ax7.plot(t * 1000, max_disp * 1e6, 'g-', linewidth=2)
    ax7.fill_between(t * 1000, 0, max_disp * 1e6, alpha=0.3, color='green')
    ax7.set_xlabel('Time (ms)', fontsize=10)
    ax7.set_ylabel('Max |u| (μm)', fontsize=10)
    ax7.set_title('Maximum Displacement Envelope', fontsize=11, fontweight='bold')
    ax7.grid(True, alpha=0.3)

    # Plot 8: RMS displacement per DOF
    ax8 = fig.add_subplot(gs[3, 2])
    rms_disp = np.sqrt(np.mean(u_time**2, axis=1))
    dof_indices = np.arange(len(rms_disp))
    ax8.bar(dof_indices, rms_disp * 1e6, alpha=0.7, edgecolor='black', linewidth=0.5)

    if contact_params and len(contact_params['contact_dofs']) > 0:
        contact_dofs = contact_params['contact_dofs']
        ax8.bar(contact_dofs, rms_disp[contact_dofs] * 1e6,
               color='red', alpha=0.7, label='Contact DOFs')
        ax8.legend(fontsize=8)

    ax8.set_xlabel('DOF Index', fontsize=10)
    ax8.set_ylabel('RMS Displacement (μm)', fontsize=10)
    ax8.set_title('RMS Displacement Distribution', fontsize=11, fontweight='bold')
    ax8.grid(True, alpha=0.3, axis='y')

    # Add overall title
    freq_hz = omega / (2 * np.pi)
    fig.suptitle(f'COPV Displacement Field Analysis: {mesh_name}\n' +
                f'Harmonic Excitation at {freq_hz:.1f} Hz',
                fontsize=14, fontweight='bold', y=0.98)

    # Save figure
    output_path = output_dir / f"copv_displacement_{mesh_name}.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n  ✓ Saved: {output_path}")

    plt.close(fig)


def main():
    """Run complete ROM nonlinear analysis for COPV meshes."""

    print("\n" + "="*70)
    print("COPV NONLINEAR ROM ANALYSIS")
    print("Craig-Bampton Reduction + Harmonic Balance AFT")
    print("="*70)

    if not ROM_AVAILABLE:
        print("\n⚠️  This analysis requires dolfinx")
        print("Install with: conda install -c conda-forge fenics-dolfinx")
        return

    # Setup output directory
    output_dir = Path("examples/output")
    output_dir.mkdir(parents=True, exist_ok=True)

    # Configuration
    n_modes = 15  # Number of CB modes per substructure

    # Analyze both configurations
    configurations = [
        {
            'file': 'bonded_cylinder_basic.msh',
            'name': 'ideal',
            'title': 'Ideal COPV',
            'has_disbond': False
        },
        {
            'file': 'bonded_cylinder_disbond.msh',
            'name': 'disbond',
            'title': 'Disbonded COPV',
            'has_disbond': True
        }
    ]

    results = {}

    for config in configurations:
        mesh_file = config['file']

        if not Path(mesh_file).exists():
            print(f"\n⚠️  Mesh file not found: {mesh_file}")
            print(f"   Run examples 04/05 first to generate meshes")
            continue

        try:
            # Step 1: Craig-Bampton Reduction
            K_r, M_r, cb = perform_rom_reduction(mesh_file, n_modes=n_modes)

            # Step 2: Nonlinear Harmonic Balance
            U_freq, u_time, f_nl_time, hb, contact_params, omega = run_nonlinear_hb_analysis(
                K_r, M_r, config['name'], has_disbond=config['has_disbond']
            )

            # Step 3: Visualize displacement fields
            visualize_displacement_snapshots(
                u_time, hb, omega, config['name'], output_dir, contact_params, f_nl_time
            )

            # Store results
            results[config['name']] = {
                'K_r': K_r,
                'M_r': M_r,
                'U_freq': U_freq,
                'u_time': u_time,
                'f_nl_time': f_nl_time,
                'omega': omega
            }

        except Exception as e:
            print(f"\n✗ Error analyzing {mesh_file}: {e}")
            import traceback
            traceback.print_exc()

    # Final summary
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)

    if results:
        print("\nGenerated visualizations:")
        for config in configurations:
            if config['name'] in results:
                print(f"  ✓ {config['title']}: copv_displacement_{config['name']}.png")

        print(f"\nAll outputs saved to: {output_dir}/")

        # Comparison
        if 'ideal' in results and 'disbond' in results:
            print("\n" + "="*70)
            print("COMPARISON: Ideal vs Disbonded COPV")
            print("="*70)

            ideal_max = np.max(np.abs(results['ideal']['u_time'])) * 1e6
            disbond_max = np.max(np.abs(results['disbond']['u_time'])) * 1e6

            print(f"\nMax Displacement:")
            print(f"  Ideal:    {ideal_max:.3f} μm")
            print(f"  Disbond:  {disbond_max:.3f} μm")
            print(f"  Ratio:    {disbond_max/ideal_max:.2f}x")

            if np.max(results['disbond']['f_nl_time']) > 0:
                max_contact = np.max(results['disbond']['f_nl_time'])
                print(f"\nContact Forces (Disbonded COPV):")
                print(f"  Max contact force: {max_contact:.3e} N")

    print("\n✅ ROM nonlinear analysis complete!")


if __name__ == "__main__":
    main()
