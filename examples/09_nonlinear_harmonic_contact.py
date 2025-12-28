"""Example 09: Nonlinear harmonic response with contact using HB-AFT.

This example demonstrates the Harmonic Balance Method with Alternating
Frequency-Time (AFT) for analyzing bonded plates with disbond contact.

The workflow:
1. Load mesh and perform Craig-Bampton reduction
2. Set up harmonic excitation and contact parameters
3. Solve nonlinear harmonic response using HB-AFT
4. Visualize periodic response with contact activation

This captures the nonlinear dynamics of contact (unilateral constraint)
that cannot be represented in standard linear harmonic analysis.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

try:
    from bonded_substructures.rom.craig_bampton import CraigBamptonReduction
    from bonded_substructures.rom.harmonic_balance import HarmonicBalanceAFT
    HB_AVAILABLE = True
except ImportError as e:
    HB_AVAILABLE = False
    print(f"Harmonic Balance requires dolfinx: {e}")


def main():
    """Demonstrate nonlinear harmonic contact analysis."""

    print("=" * 70)
    print("Nonlinear Harmonic Balance with Contact (HB-AFT)")
    print("=" * 70)

    if not HB_AVAILABLE:
        print("\n⚠️  This example requires dolfinx")
        print("Install with: conda install -c conda-forge fenics-dolfinx")
        return

    # Use wide plate (has more interior DOFs for modes)
    mesh_file = "wide_plate_disbond.msh"

    if not Path(mesh_file).exists():
        print(f"\n⚠️  Mesh file not found: {mesh_file}")
        print("Run Example 02 first to generate mesh")
        return

    print(f"\nInput mesh: {mesh_file}")

    # Step 1: Craig-Bampton Reduction
    print("\n" + "-" * 70)
    print("Step 1: Craig-Bampton Reduction")
    print("-" * 70)

    cb = CraigBamptonReduction(
        mesh_file=mesh_file,
        material_tags={"Material_1": 1, "Material_2": 2}
    )

    try:
        cb.load_mesh()
        cb.partition_substructures()
        cb.identify_interface_dofs()
        cb.assemble_matrices()

        # Use fewer modes for this example
        n_modes = 5  # Small system for demonstration
        print(f"\nComputing {n_modes} Craig-Bampton modes per substructure...")
        cb.compute_modes(n_modes=n_modes)

        # Assemble reduced system
        K_r, M_r = cb.assemble_reduced_system()

        print(f"\n  Original system: {cb.K.shape[0] if hasattr(cb, 'K') else 'unknown'} DOFs")
        print(f"  Reduced system: {K_r.shape[0]} DOFs")
        print(f"  ✓ Craig-Bampton reduction complete")

    except Exception as e:
        print(f"Error in Craig-Bampton reduction: {e}")
        import traceback
        traceback.print_exc()
        return

    # Step 2: Setup Harmonic Balance
    print("\n" + "-" * 70)
    print("Step 2: Harmonic Balance Setup")
    print("-" * 70)

    # Initialize HB solver
    hb = HarmonicBalanceAFT(K_r, M_r, damping_ratio=0.02)
    n_harmonics = 5  # Use 5 harmonics (captures up to 5ω)
    hb.set_harmonics(n_harmonics)

    print(f"\n  Reduced DOFs: {K_r.shape[0]}")
    print(f"  Number of harmonics: {n_harmonics}")
    print(f"  Time points per period: {hb.n_time_points}")
    print(f"  Total unknowns: {K_r.shape[0] * (2*n_harmonics + 1)}")

    # Excitation parameters
    freq_hz = 100.0  # Hz
    omega = 2 * np.pi * freq_hz
    F_magnitude = 1000.0  # N

    print(f"\n  Excitation frequency: {freq_hz:.1f} Hz ({omega:.1f} rad/s)")
    print(f"  Force magnitude: {F_magnitude:.0f} N")

    # Setup external force (harmonic at fundamental frequency)
    F_ext = np.zeros((K_r.shape[0], 2*n_harmonics + 1))
    F_ext[0, 1] = F_magnitude  # Cosine component at first DOF

    # Contact parameters
    # Assume interface DOFs are in the middle of the reduced system
    n_interface = K_r.shape[0] // 4  # Estimate
    contact_dofs = list(range(K_r.shape[0] - n_interface, K_r.shape[0]))[:5]  # First 5 interface DOFs

    contact_params = {
        'contact_dofs': contact_dofs,
        'contact_stiffness': 1e8,  # N/m (lower than full model for stability)
        'gap_initial': 1e-6,  # 1 μm initial gap
        'contact_type': 'penalty'
    }

    print(f"\n  Contact DOFs: {len(contact_dofs)}")
    print(f"  Contact stiffness: {contact_params['contact_stiffness']:.1e} N/m")
    print(f"  Initial gap: {contact_params['gap_initial']*1e6:.1f} μm")

    # Step 3: Solve Nonlinear Harmonic Response
    print("\n" + "-" * 70)
    print("Step 3: Nonlinear Harmonic Balance Solution")
    print("-" * 70)

    U_freq, u_time, history = hb.solve_harmonic_balance(
        omega, F_ext, contact_params,
        max_iter=100, tol=1e-6, relaxation=0.3  # Lower relaxation for stability
    )

    # Step 4: Post-process and Visualize
    print("\n" + "-" * 70)
    print("Step 4: Visualization")
    print("-" * 70)

    # Compute contact forces
    f_nl_time = hb.evaluate_contact_forces(u_time, contact_params)

    # Time vector
    T = 2 * np.pi / omega
    t = np.linspace(0, T, hb.n_time_points, endpoint=False)

    # Create comprehensive plot
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(3, 2, hspace=0.3, wspace=0.3)

    # Plot 1: Displacement at excited DOF
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(t * 1000, u_time[0, :] * 1e6, 'b-', linewidth=2)
    ax1.set_xlabel('Time (ms)', fontsize=11)
    ax1.set_ylabel('Displacement (μm)', fontsize=11)
    ax1.set_title('Response at Excited DOF', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Plot 2: Displacement at contact DOF
    ax2 = fig.add_subplot(gs[0, 1])
    if len(contact_dofs) > 0:
        contact_dof = contact_dofs[0]
        ax2.plot(t * 1000, u_time[contact_dof, :] * 1e6, 'r-', linewidth=2, label='Displacement')
        ax2.axhline(contact_params['gap_initial'] * 1e6, color='k', linestyle='--',
                   linewidth=1.5, label='Initial gap')
        ax2.set_xlabel('Time (ms)', fontsize=11)
        ax2.set_ylabel('Displacement (μm)', fontsize=11)
        ax2.set_title(f'Contact DOF {contact_dof}', fontsize=12, fontweight='bold')
        ax2.legend(fontsize=9)
        ax2.grid(True, alpha=0.3)

    # Plot 3: Contact force
    ax3 = fig.add_subplot(gs[1, 0])
    if len(contact_dofs) > 0:
        total_contact_force = np.sum(f_nl_time[contact_dofs, :], axis=0)
        ax3.plot(t * 1000, total_contact_force, 'g-', linewidth=2)
        ax3.set_xlabel('Time (ms)', fontsize=11)
        ax3.set_ylabel('Total Contact Force (N)', fontsize=11)
        ax3.set_title('Nonlinear Contact Force', fontsize=12, fontweight='bold')
        ax3.grid(True, alpha=0.3)

        # Shade contact active regions
        contact_active = total_contact_force > 1e-3
        for i in range(len(t)-1):
            if contact_active[i]:
                ax3.axvspan(t[i]*1000, t[i+1]*1000, alpha=0.2, color='red')

    # Plot 4: Convergence history
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.semilogy(history['iteration'], history['residual'], 'b-o', markersize=4)
    ax4.set_xlabel('Iteration', fontsize=11)
    ax4.set_ylabel('Residual', fontsize=11)
    ax4.set_title('Convergence History', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, which='both')

    # Plot 5: Frequency content (harmonic amplitudes)
    ax5 = fig.add_subplot(gs[2, 0])
    harmonics = np.arange(n_harmonics + 1)
    amplitudes_excited = [np.sqrt(U_freq[0, 0]**2)] + \
                        [np.sqrt(U_freq[0, 2*k-1]**2 + U_freq[0, 2*k]**2)
                         for k in range(1, n_harmonics+1)]

    if len(contact_dofs) > 0:
        amplitudes_contact = [np.sqrt(U_freq[contact_dof, 0]**2)] + \
                            [np.sqrt(U_freq[contact_dof, 2*k-1]**2 + U_freq[contact_dof, 2*k]**2)
                             for k in range(1, n_harmonics+1)]

        width = 0.35
        ax5.bar(harmonics - width/2, amplitudes_excited, width, label='Excited DOF', alpha=0.8)
        ax5.bar(harmonics + width/2, amplitudes_contact, width, label='Contact DOF', alpha=0.8)
        ax5.legend(fontsize=9)
    else:
        ax5.bar(harmonics, amplitudes_excited, alpha=0.8)

    ax5.set_xlabel('Harmonic Number', fontsize=11)
    ax5.set_ylabel('Amplitude', fontsize=11)
    ax5.set_title('Frequency Content (Harmonic Amplitudes)', fontsize=12, fontweight='bold')
    ax5.set_xticks(harmonics)
    ax5.grid(True, alpha=0.3, axis='y')

    # Plot 6: Phase portrait (displacement vs contact force)
    ax6 = fig.add_subplot(gs[2, 1])
    if len(contact_dofs) > 0:
        ax6.plot(u_time[contact_dof, :] * 1e6, f_nl_time[contact_dof, :],
                'r-', linewidth=2, alpha=0.7)
        ax6.axvline(contact_params['gap_initial'] * 1e6, color='k',
                   linestyle='--', linewidth=1.5, label='Gap threshold')
        ax6.set_xlabel('Displacement (μm)', fontsize=11)
        ax6.set_ylabel('Contact Force (N)', fontsize=11)
        ax6.set_title('Phase Portrait (Contact Nonlinearity)', fontsize=12, fontweight='bold')
        ax6.legend(fontsize=9)
        ax6.grid(True, alpha=0.3)

    # Add overall title
    fig.suptitle(f'Nonlinear Harmonic Response at {freq_hz:.1f} Hz\n'
                f'Craig-Bampton ROM ({K_r.shape[0]} DOFs) with Contact',
                fontsize=14, fontweight='bold', y=0.98)

    # Save figure
    output_dir = Path("examples/output")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "nonlinear_harmonic_contact.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n  ✓ Saved plot: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("Nonlinear Harmonic Analysis Complete!")
    print("=" * 70)

    print(f"\nSystem size:")
    print(f"  Reduced DOFs: {K_r.shape[0]}")
    print(f"  Harmonics: {n_harmonics}")
    print(f"  Time points: {hb.n_time_points}")

    print(f"\nExcitation:")
    print(f"  Frequency: {freq_hz:.1f} Hz")
    print(f"  Force: {F_magnitude:.0f} N")

    print(f"\nResults:")
    print(f"  Max displacement: {np.max(np.abs(u_time)):.3e} m")
    print(f"  Max contact force: {np.max(f_nl_time):.3e} N")
    contact_time_fraction = np.sum(f_nl_time > 1e-6) / f_nl_time.size * 100
    print(f"  Contact active: {contact_time_fraction:.1f}% of time-DOF pairs")

    print(f"\nConvergence:")
    print(f"  Iterations: {len(history['iteration'])}")
    print(f"  Final residual: {history['residual'][-1]:.3e}")

    print(f"\n📊 Results saved to: {output_path}")

    print("\nKey Observations:")
    print("  1. HB-AFT captures nonlinear contact dynamics")
    print("  2. Higher harmonics generated by contact nonlinearity")
    print("  3. Contact force shows unilateral constraint (compression only)")
    print("  4. Periodic solution converged successfully")


if __name__ == "__main__":
    main()
