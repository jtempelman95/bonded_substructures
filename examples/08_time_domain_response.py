"""Example 08: Time-domain response using Craig-Bampton ROM.

This example demonstrates transient dynamic analysis using the
Craig-Bampton reduced order model with Newmark-beta time integration.

COORDINATE SYSTEM:
- x-y plane: In-plane dimensions (width × length)
- z-direction: Through-thickness (material stacking)
- Bond interface: x-y plane at z = t1

Loading scenarios:
1. Impulse load (impact)
2. Step load (sudden force application)
3. Harmonic load (sinusoidal excitation)

The ROM enables efficient time-domain simulation of large structural systems.

Mesh files:
- Uses wide_plate_disbond.msh from Example 07
- Falls back to bonded_plate_disbond.msh from Example 02 if not available
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


def newmark_beta(M, C, K, F_func, dt, t_end, u0=None, v0=None,
                 beta=0.25, gamma=0.5):
    """Newmark-beta time integration for second-order systems.

    Solves: M*ü + C*u̇ + K*u = F(t)

    Args:
        M: Mass matrix (n x n)
        C: Damping matrix (n x n)
        K: Stiffness matrix (n x n)
        F_func: Function F(t) returning force vector
        dt: Time step (s)
        t_end: End time (s)
        u0: Initial displacement (default: zeros)
        v0: Initial velocity (default: zeros)
        beta: Newmark beta parameter (default: 0.25 for constant acceleration)
        gamma: Newmark gamma parameter (default: 0.5 for constant acceleration)

    Returns:
        t: Time array
        u: Displacement history (n_dof x n_time)
        v: Velocity history (n_dof x n_time)
        a: Acceleration history (n_dof x n_time)
    """
    n_dof = M.shape[0]
    n_steps = int(t_end / dt) + 1
    t = np.linspace(0, t_end, n_steps)

    # Initial conditions
    if u0 is None:
        u0 = np.zeros(n_dof)
    if v0 is None:
        v0 = np.zeros(n_dof)

    # Initialize arrays
    u = np.zeros((n_dof, n_steps))
    v = np.zeros((n_dof, n_steps))
    a = np.zeros((n_dof, n_steps))

    u[:, 0] = u0
    v[:, 0] = v0

    # Initial acceleration: a0 = M^(-1) * (F0 - C*v0 - K*u0)
    F0 = F_func(t[0])
    a[:, 0] = np.linalg.solve(M, F0 - C @ v0 - K @ u0)

    # Newmark integration constants
    a0 = 1.0 / (beta * dt**2)
    a1 = gamma / (beta * dt)
    a2 = 1.0 / (beta * dt)
    a3 = 1.0 / (2 * beta) - 1.0
    a4 = gamma / beta - 1.0
    a5 = dt / 2 * (gamma / beta - 2.0)
    a6 = dt * (1.0 - gamma)
    a7 = gamma * dt

    # Effective stiffness matrix
    K_eff = K + a0 * M + a1 * C

    print(f"\n  Time integration:")
    print(f"    Time step: {dt:.4f} s")
    print(f"    Total time: {t_end:.2f} s")
    print(f"    Number of steps: {n_steps}")
    print(f"    Newmark parameters: beta={beta}, gamma={gamma}")

    # Time stepping
    for i in range(n_steps - 1):
        # Effective force
        F_eff = (F_func(t[i+1]) +
                 M @ (a0 * u[:, i] + a2 * v[:, i] + a3 * a[:, i]) +
                 C @ (a1 * u[:, i] + a4 * v[:, i] + a5 * a[:, i]))

        # Solve for displacement at next time step
        u[:, i+1] = np.linalg.solve(K_eff, F_eff)

        # Update acceleration and velocity
        a[:, i+1] = a0 * (u[:, i+1] - u[:, i]) - a2 * v[:, i] - a3 * a[:, i]
        v[:, i+1] = v[:, i] + a6 * a[:, i] + a7 * a[:, i+1]

        if (i + 1) % (n_steps // 10) == 0:
            print(f"      Progress: {(i+1)/n_steps*100:.0f}%")

    print(f"    ✓ Time integration complete")

    return t, u, v, a


def main():
    """Demonstrate time-domain response with Craig-Bampton ROM."""

    print("=" * 70)
    print("Craig-Bampton ROM - Time-Domain Response")
    print("=" * 70)

    if not CB_AVAILABLE:
        print("\n⚠️  Craig-Bampton reduction requires dolfinx")
        print("Install with: conda install -c conda-forge fenics-dolfinx")
        return

    # Check for wide plate mesh first, fall back to basic plate if not available
    mesh_file = "wide_plate_disbond.msh"
    if not Path(mesh_file).exists():
        print(f"\n⚠️  Wide plate mesh not found, using basic plate mesh")
        mesh_file = "bonded_plate_disbond.msh"

    if not Path(mesh_file).exists():
        print(f"\n⚠️  Mesh file not found: {mesh_file}")
        print("Run Example 02 or 07 first to generate mesh")
        return

    print(f"\nInput mesh: {mesh_file}")

    # Craig-Bampton reduction
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

    # Time-domain simulation parameters
    print("\n" + "-" * 70)
    print("Step 2: Time-Domain Simulation")
    print("-" * 70)

    # Simulation parameters
    dt = 0.0001  # Time step (s) - 0.1 ms
    t_end = 0.05  # Total time (s) - 50 ms
    damping_ratio = 0.02  # 2% damping

    # Rayleigh damping: C = alpha*M + beta*K
    # Estimate modal frequencies using generalized eigenvalue problem
    from scipy.linalg import eigh
    eigenvals, _ = eigh(K_r, M_r)
    eigenvals = eigenvals[eigenvals > 0]  # Keep only positive eigenvalues
    omega_min = np.sqrt(eigenvals.min())
    omega_max = np.sqrt(eigenvals.max())

    alpha = 2 * damping_ratio * omega_min * omega_max / (omega_min + omega_max)
    beta = 2 * damping_ratio / (omega_min + omega_max)

    C_r = alpha * M_r + beta * K_r

    print(f"\n  System properties:")
    print(f"    Reduced DOFs: {K_r.shape[0]}")
    print(f"    Modal frequency range: {omega_min/(2*np.pi):.1f} - {omega_max/(2*np.pi):.1f} Hz")
    print(f"    Damping: {damping_ratio*100:.1f}% (Rayleigh)")

    # Define loading scenarios
    print(f"\n  Loading scenario: Impulse load (impact)")

    # Impulse load (short duration impact)
    F_magnitude = 10000.0  # N
    impulse_duration = 0.001  # 1 ms
    force_dof = 0  # Apply to first DOF

    def impulse_force(t):
        """Impulse force (triangular pulse)."""
        F = np.zeros(K_r.shape[0])
        if t <= impulse_duration:
            # Triangular pulse
            F[force_dof] = F_magnitude * (1 - t / impulse_duration)
        return F

    # Time integration
    print(f"\n  Impulse parameters:")
    print(f"    Magnitude: {F_magnitude:.0f} N")
    print(f"    Duration: {impulse_duration*1000:.2f} ms")
    print(f"    Location: DOF {force_dof}")

    t, u, v, a = newmark_beta(M_r, C_r, K_r, impulse_force, dt, t_end)

    # Extract response at observation DOF
    obs_dof = force_dof
    displacement = u[obs_dof, :]
    velocity = v[obs_dof, :]
    acceleration = a[obs_dof, :]

    # Statistics
    print(f"\n  Response statistics:")
    print(f"    Max displacement: {np.max(np.abs(displacement)):.2e} m")
    print(f"    Max velocity: {np.max(np.abs(velocity)):.2e} m/s")
    print(f"    Max acceleration: {np.max(np.abs(acceleration)):.2e} m/s²")

    # Visualization
    print("\n" + "-" * 70)
    print("Step 3: Visualization")
    print("-" * 70)

    fig, axes = plt.subplots(4, 1, figsize=(12, 14))

    # Plot 1: Force history
    ax1 = axes[0]
    force_history = np.array([impulse_force(ti)[force_dof] for ti in t])
    ax1.plot(t * 1000, force_history / 1000, 'r-', linewidth=2)
    ax1.set_ylabel('Force (kN)', fontsize=12)
    ax1.set_title('Input: Impulse Load', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, t_end * 1000])

    # Plot 2: Displacement response
    ax2 = axes[1]
    ax2.plot(t * 1000, displacement * 1000, 'b-', linewidth=1.5)
    ax2.set_ylabel('Displacement (mm)', fontsize=12)
    ax2.set_title('Response: Displacement', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, t_end * 1000])

    # Plot 3: Velocity response
    ax3 = axes[2]
    ax3.plot(t * 1000, velocity * 1000, 'g-', linewidth=1.5)
    ax3.set_ylabel('Velocity (mm/s)', fontsize=12)
    ax3.set_title('Response: Velocity', fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([0, t_end * 1000])

    # Plot 4: Acceleration response
    ax4 = axes[3]
    ax4.plot(t * 1000, acceleration, 'm-', linewidth=1.5)
    ax4.set_xlabel('Time (ms)', fontsize=12)
    ax4.set_ylabel('Acceleration (m/s²)', fontsize=12)
    ax4.set_title('Response: Acceleration', fontsize=14, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.set_xlim([0, t_end * 1000])

    plt.tight_layout()

    # Save figure
    output_dir = Path("examples/output")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "time_domain_response_rom.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n  ✓ Saved plot: {output_path}")

    # Summary
    print("\n" + "=" * 70)
    print("Time-Domain Analysis Complete!")
    print("=" * 70)
    print(f"\nSimulation parameters:")
    print(f"  Time step: {dt*1000:.2f} ms")
    print(f"  Total time: {t_end*1000:.1f} ms")
    print(f"  Number of steps: {len(t)}")
    print(f"  ROM size: {K_r.shape[0]} DOFs (reduction: {690/K_r.shape[0]:.1f}x)")

    print(f"\n📊 Results saved to: {output_path}")


if __name__ == "__main__":
    main()
