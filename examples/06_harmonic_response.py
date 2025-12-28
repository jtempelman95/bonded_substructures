"""Example 06: Harmonic response using Craig-Bampton ROM.

This example demonstrates the complete Craig-Bampton workflow:
1. Load mesh and partition into substructures
2. Identify interior and interface DOFs
3. Assemble mass and stiffness matrices
4. Compute constraint modes and fixed-interface normal modes
5. Assemble reduced system (50 modes per substructure)
6. Apply contact constraints at interface
7. Simulate harmonic response to sinusoidal loading

The reduced system dramatically reduces computational cost while
preserving accuracy for dynamic response.
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

try:
    from bonded_substructures.rom.craig_bampton import CraigBamptonReduction
    CB_AVAILABLE = True
except ImportError as e:
    CB_AVAILABLE = False
    print(f"Craig-Bampton functionality not available: {e}")
    print("Install dolfinx to use this example")


def harmonic_response(K, M, F_func, omega_range, damping_ratio=0.02):
    """Compute frequency response function.

    Args:
        K: Stiffness matrix (n x n)
        M: Mass matrix (n x n)
        F_func: Function that returns force vector for each frequency
        omega_range: Array of angular frequencies (rad/s)
        damping_ratio: Modal damping ratio (default 2%)

    Returns:
        response: Displacement response array (n_dof x n_freq)
    """
    n_dof = K.shape[0]
    n_freq = len(omega_range)
    response = np.zeros((n_dof, n_freq), dtype=complex)

    # Proportional damping: C = alpha*M + beta*K
    # For modal damping: zeta = (alpha/(2*omega) + beta*omega/2)
    # Simplification: use Rayleigh damping with beta = 0
    omega_ref = omega_range[len(omega_range)//2]  # Reference frequency
    alpha = 2 * damping_ratio * omega_ref
    beta = 0.0

    C = alpha * M + beta * K

    print(f"\n  Solving frequency response...")
    print(f"    Frequency range: {omega_range[0]/(2*np.pi):.1f} - {omega_range[-1]/(2*np.pi):.1f} Hz")
    print(f"    Number of frequencies: {n_freq}")
    print(f"    Damping ratio: {damping_ratio*100:.1f}%")

    for i, omega in enumerate(omega_range):
        # Dynamic stiffness matrix: K_dyn = K - omega^2 * M + 1j * omega * C
        K_dyn = K - omega**2 * M + 1j * omega * C

        # Force vector at this frequency
        F = F_func(omega)

        # Solve: K_dyn * u = F
        u = np.linalg.solve(K_dyn, F)
        response[:, i] = u

        if (i + 1) % 20 == 0:
            print(f"      Progress: {i+1}/{n_freq} frequencies")

    print(f"    ✓ Frequency response complete")

    return response


def main():
    """Demonstrate Craig-Bampton ROM with harmonic response."""

    print("=" * 70)
    print("Craig-Bampton ROM - Harmonic Response Simulation")
    print("=" * 70)

    if not CB_AVAILABLE:
        print("\n⚠️  Craig-Bampton reduction requires dolfinx")
        print("Install with: conda install -c conda-forge fenics-dolfinx")
        return

    # Use conformal mesh
    mesh_file = "rectangle_disbond.msh"

    if not Path(mesh_file).exists():
        print(f"\n⚠️  Mesh file not found: {mesh_file}")
        print("Run Example 02 first to generate disbond mesh")
        return

    print(f"\nInput mesh: {mesh_file}")

    # Initialize Craig-Bampton reduction
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

        # Compute modes (50 per substructure)
        n_modes = 50
        cb.compute_modes(n_modes=n_modes)

        # Assemble reduced system
        K_r, M_r = cb.assemble_reduced_system()

        # Apply contact constraints
        cb.apply_contact_constraints(contact_stiffness=1e10)

    except Exception as e:
        print(f"Error in Craig-Bampton reduction: {e}")
        import traceback
        traceback.print_exc()
        return

    # Harmonic response analysis
    print("\n" + "-" * 70)
    print("Step 2: Harmonic Response Analysis")
    print("-" * 70)

    # Frequency range for analysis
    f_min = 1.0  # Hz
    f_max = 500.0  # Hz
    n_freq = 200
    frequencies = np.linspace(f_min, f_max, n_freq)
    omega_range = 2 * np.pi * frequencies

    # Define harmonic force (point load at first DOF)
    F_amplitude = 1000.0  # N
    force_dof = 0  # Apply force to first DOF

    def force_vector(omega):
        """Force vector as function of frequency."""
        F = np.zeros(K_r.shape[0])
        F[force_dof] = F_amplitude
        return F

    # Compute frequency response
    response = harmonic_response(K_r, M_r, force_vector, omega_range, damping_ratio=0.02)

    # Extract response at observation DOF
    obs_dof = force_dof
    displacement = np.abs(response[obs_dof, :])

    # Find resonance peaks
    peaks_idx = []
    for i in range(1, len(displacement) - 1):
        if displacement[i] > displacement[i-1] and displacement[i] > displacement[i+1]:
            if displacement[i] > 0.1 * np.max(displacement):  # Significant peaks only
                peaks_idx.append(i)

    print(f"\n  Response Statistics:")
    print(f"    Maximum displacement: {np.max(displacement):.2e} m")
    print(f"    Resonance peaks detected: {len(peaks_idx)}")
    if peaks_idx:
        print(f"    Resonance frequencies:")
        for idx in peaks_idx[:5]:  # Show first 5 peaks
            print(f"      {frequencies[idx]:.1f} Hz (displacement: {displacement[idx]:.2e} m)")

    # Visualization
    print("\n" + "-" * 70)
    print("Step 3: Visualization")
    print("-" * 70)

    fig, axes = plt.subplots(2, 1, figsize=(12, 10))

    # Plot 1: Frequency Response Function
    ax1 = axes[0]
    ax1.semilogy(frequencies, displacement, 'b-', linewidth=1.5, label='Displacement amplitude')
    ax1.scatter(frequencies[peaks_idx], displacement[peaks_idx],
                c='r', s=100, zorder=5, label='Resonance peaks')
    ax1.set_xlabel('Frequency (Hz)', fontsize=12)
    ax1.set_ylabel('Displacement Amplitude (m)', fontsize=12)
    ax1.set_title(f'Frequency Response Function (ROM with {n_modes} modes/substructure)', fontsize=14)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.set_xlim([f_min, f_max])

    # Plot 2: Phase Response
    ax2 = axes[1]
    phase = np.angle(response[obs_dof, :], deg=True)
    ax2.plot(frequencies, phase, 'g-', linewidth=1.5)
    ax2.set_xlabel('Frequency (Hz)', fontsize=12)
    ax2.set_ylabel('Phase (degrees)', fontsize=12)
    ax2.set_title('Phase Response', fontsize=14)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([f_min, f_max])
    ax2.set_ylim([-180, 180])

    plt.tight_layout()

    # Save figure
    output_dir = Path("examples/output")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "harmonic_response_rom.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n  ✓ Saved plot: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("Harmonic Response Analysis Complete!")
    print("=" * 70)
    print(f"\nReduced system size: {K_r.shape[0]} DOFs")
    print(f"  Original system: 690 DOFs")
    print(f"  Reduction factor: {690/K_r.shape[0]:.1f}x")
    print(f"\nModal content:")
    for name, sub in cb.substructures.items():
        if sub.frequencies is not None:
            print(f"  {name}: {len(sub.frequencies)} modes, "
                  f"range {sub.frequencies[0]:.1f} - {sub.frequencies[-1]:.1f} Hz")

    print(f"\n📊 Results saved to: {output_path}")


if __name__ == "__main__":
    main()
