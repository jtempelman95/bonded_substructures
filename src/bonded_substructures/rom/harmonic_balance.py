"""Harmonic Balance Method with Alternating Frequency-Time (HB-AFT).

This module implements nonlinear harmonic analysis for structures with
contact nonlinearities using the Harmonic Balance Method with AFT.

Theory:
-------
For a nonlinear system with contact:
    M*ü + C*u̇ + K*u + f_nl(u, u̇) = F*cos(ωt)

where f_nl represents nonlinear contact forces.

Harmonic Balance assumes a periodic solution:
    u(t) = Σ [A_k*cos(kωt) + B_k*sin(kωt)]
         k=0..N_h

AFT Strategy:
1. Represent displacement as Fourier series (frequency domain)
2. Use FFT to convert to time domain
3. Evaluate nonlinear forces in time domain
4. Use IFFT to return to frequency domain
5. Iterate until convergence

This is extremely efficient for contact because:
- Linear dynamics in frequency domain (exact)
- Nonlinear contact in time domain (accurate)
- FFT makes conversion O(N log N)

References:
-----------
[1] Cameron & Griffin, "An alternating frequency/time domain method for
    calculating the steady-state response of nonlinear dynamic systems,"
    J. Applied Mechanics, 1989.
[2] Krack & Groß, "Harmonic Balance for Nonlinear Vibration Problems," 2019.
"""

import numpy as np
from scipy.linalg import solve, eigh
from scipy.fft import fft, ifft, fftfreq
import warnings


class HarmonicBalanceAFT:
    """Harmonic Balance solver with Alternating Frequency-Time method.

    Solves nonlinear periodic response of reduced-order models with
    contact nonlinearities.

    Parameters
    ----------
    K : ndarray
        Stiffness matrix (n_dof × n_dof)
    M : ndarray
        Mass matrix (n_dof × n_dof)
    C : ndarray, optional
        Damping matrix (n_dof × n_dof). If None, uses Rayleigh damping.
    damping_ratio : float, optional
        Modal damping ratio (default: 0.02 = 2%)

    Attributes
    ----------
    n_dof : int
        Number of degrees of freedom
    n_harmonics : int
        Number of harmonics in Fourier series
    n_time_points : int
        Number of time points (2*n_harmonics + 1)
    """

    def __init__(self, K, M, C=None, damping_ratio=0.02):
        """Initialize Harmonic Balance solver."""
        self.K = K
        self.M = M
        self.n_dof = K.shape[0]
        self.damping_ratio = damping_ratio

        # Compute damping matrix if not provided
        if C is None:
            # Rayleigh damping: C = α*M + β*K
            eigenvals, _ = eigh(K, M)
            eigenvals = eigenvals[eigenvals > 0]
            omega_min = np.sqrt(eigenvals.min())
            omega_max = np.sqrt(eigenvals.max())

            alpha = 2 * damping_ratio * omega_min * omega_max / (omega_min + omega_max)
            beta = 2 * damping_ratio / (omega_min + omega_max)

            self.C = alpha * M + beta * K
        else:
            self.C = C

    def set_harmonics(self, n_harmonics):
        """Set number of harmonics for Fourier series.

        Parameters
        ----------
        n_harmonics : int
            Number of harmonics (excluding DC component)
            Total unknowns = (2*n_harmonics + 1) * n_dof
        """
        self.n_harmonics = n_harmonics
        self.n_time_points = 2 * n_harmonics + 1

    def frequency_to_time(self, U_freq, omega):
        """Convert frequency-domain solution to time domain using IFFT.

        Parameters
        ----------
        U_freq : ndarray
            Frequency domain solution (n_dof × (2*n_h+1))
            Format: [A_0, A_1, B_1, A_2, B_2, ..., A_n, B_n]
        omega : float
            Excitation frequency (rad/s)

        Returns
        -------
        u_time : ndarray
            Time-domain displacement (n_dof × n_time_points)
        t : ndarray
            Time points over one period
        """
        T = 2 * np.pi / omega
        t = np.linspace(0, T, self.n_time_points, endpoint=False)

        n_dof = U_freq.shape[0]
        u_time = np.zeros((n_dof, self.n_time_points))

        # Reconstruct time series from Fourier coefficients
        # U_freq format: [A_0, A_1, B_1, A_2, B_2, ...]
        for i in range(n_dof):
            # DC component
            u_time[i, :] = U_freq[i, 0]

            # Harmonic components
            for k in range(1, self.n_harmonics + 1):
                A_k = U_freq[i, 2*k - 1]  # Cosine coefficient
                B_k = U_freq[i, 2*k]      # Sine coefficient
                u_time[i, :] += A_k * np.cos(k * omega * t) + B_k * np.sin(k * omega * t)

        return u_time, t

    def time_to_frequency(self, f_nl_time, omega):
        """Convert time-domain nonlinear forces to frequency domain using FFT.

        Parameters
        ----------
        f_nl_time : ndarray
            Nonlinear forces in time domain (n_dof × n_time_points)
        omega : float
            Excitation frequency (rad/s)

        Returns
        -------
        F_nl_freq : ndarray
            Frequency domain forces (n_dof × (2*n_h+1))
            Format: [A_0, A_1, B_1, A_2, B_2, ..., A_n, B_n]
        """
        n_dof = f_nl_time.shape[0]
        F_nl_freq = np.zeros((n_dof, 2 * self.n_harmonics + 1))

        T = 2 * np.pi / omega

        for i in range(n_dof):
            # Use FFT to get Fourier coefficients
            f_fft = fft(f_nl_time[i, :])

            # Extract DC component
            F_nl_freq[i, 0] = np.real(f_fft[0]) / self.n_time_points

            # Extract harmonic components
            for k in range(1, self.n_harmonics + 1):
                # Cosine coefficient (real part)
                F_nl_freq[i, 2*k - 1] = 2 * np.real(f_fft[k]) / self.n_time_points
                # Sine coefficient (imaginary part)
                F_nl_freq[i, 2*k] = -2 * np.imag(f_fft[k]) / self.n_time_points

        return F_nl_freq

    def assemble_frequency_matrix(self, omega):
        """Assemble frequency-domain dynamic stiffness matrix.

        For harmonic k at frequency ω:
            Z_k = K - (kω)²*M + i*(kω)*C

        Parameters
        ----------
        omega : float
            Fundamental frequency (rad/s)

        Returns
        -------
        Z : ndarray
            Block-diagonal frequency matrix ((2*n_h+1)*n_dof × (2*n_h+1)*n_dof)
        """
        n_total = (2 * self.n_harmonics + 1) * self.n_dof
        Z = np.zeros((n_total, n_total))

        # DC component (k=0): just stiffness
        Z[:self.n_dof, :self.n_dof] = self.K

        # Harmonic components
        for k in range(1, self.n_harmonics + 1):
            omega_k = k * omega

            # Real part: K - ω²M
            K_real = self.K - omega_k**2 * self.M
            # Imaginary part: ωC
            C_imag = omega_k * self.C

            # Indices for this harmonic
            idx_cos = slice((2*k-1)*self.n_dof, 2*k*self.n_dof)
            idx_sin = slice(2*k*self.n_dof, (2*k+1)*self.n_dof)

            # Block matrix for complex impedance
            # [Z_real  -Z_imag] [A_k]   [F_cos_k]
            # [Z_imag   Z_real] [B_k] = [F_sin_k]
            Z[idx_cos, idx_cos] = K_real
            Z[idx_cos, idx_sin] = -C_imag
            Z[idx_sin, idx_cos] = C_imag
            Z[idx_sin, idx_sin] = K_real

        return Z

    def evaluate_contact_forces(self, u_time, contact_params):
        """Evaluate nonlinear contact forces in time domain.

        Contact model (unilateral constraint):
            f_nl = k_contact * max(0, -gap)     (penalty, no penetration)
                 = 0                             (separation allowed)

        Parameters
        ----------
        u_time : ndarray
            Displacement in time domain (n_dof × n_time_points)
        contact_params : dict
            Contact parameters:
                'contact_dofs': array of DOFs with contact
                'contact_stiffness': penalty stiffness (N/m)
                'gap_initial': initial gap at each contact DOF
                'contact_type': 'penalty' or 'complementarity'

        Returns
        -------
        f_nl_time : ndarray
            Nonlinear contact forces (n_dof × n_time_points)
        """
        f_nl_time = np.zeros_like(u_time)

        contact_dofs = contact_params.get('contact_dofs', [])
        k_contact = contact_params.get('contact_stiffness', 1e10)
        gap_initial = contact_params.get('gap_initial', 0.0)
        contact_type = contact_params.get('contact_type', 'penalty')

        if len(contact_dofs) == 0:
            return f_nl_time

        # Evaluate contact at each time point
        for i, dof in enumerate(contact_dofs):
            for t_idx in range(self.n_time_points):
                # Current gap (negative = penetration)
                if isinstance(gap_initial, (int, float)):
                    gap = gap_initial - u_time[dof, t_idx]
                else:
                    gap = gap_initial[i] - u_time[dof, t_idx]

                if contact_type == 'penalty':
                    # Penalty method: restoring force when gap < 0
                    if gap < 0:
                        f_nl_time[dof, t_idx] = k_contact * abs(gap)
                    # else: no force (separation allowed)

                elif contact_type == 'complementarity':
                    # More accurate: f ≥ 0, gap ≥ 0, f*gap = 0
                    # This requires iterative solution
                    if gap < 0:
                        f_nl_time[dof, t_idx] = k_contact * abs(gap)

        return f_nl_time

    def solve_harmonic_balance(self, omega, F_ext, contact_params,
                               max_iter=50, tol=1e-6, relaxation=0.5):
        """Solve nonlinear harmonic response using AFT method.

        Algorithm (AFT iteration):
        1. Initialize: U_freq = linear solution
        2. Convert to time: u_time = IFFT(U_freq)
        3. Evaluate contact: f_nl_time = contact(u_time)
        4. Convert to frequency: F_nl_freq = FFT(f_nl_time)
        5. Update: Z*U_freq = F_ext - F_nl_freq
        6. Check convergence, repeat

        Parameters
        ----------
        omega : float
            Excitation frequency (rad/s)
        F_ext : ndarray
            External force vector in frequency domain (n_dof × (2*n_h+1))
        contact_params : dict
            Contact parameters (see evaluate_contact_forces)
        max_iter : int, optional
            Maximum iterations (default: 50)
        tol : float, optional
            Convergence tolerance (default: 1e-6)
        relaxation : float, optional
            Relaxation factor for stability (0 < r ≤ 1, default: 0.5)

        Returns
        -------
        U_freq : ndarray
            Converged frequency-domain solution (n_dof × (2*n_h+1))
        u_time : ndarray
            Time-domain solution (n_dof × n_time_points)
        history : dict
            Convergence history
        """
        print(f"\n{'='*70}")
        print("Harmonic Balance AFT Solver")
        print(f"{'='*70}")
        print(f"  Frequency: {omega/(2*np.pi):.2f} Hz ({omega:.2f} rad/s)")
        print(f"  DOFs: {self.n_dof}")
        print(f"  Harmonics: {self.n_harmonics}")
        print(f"  Time points: {self.n_time_points}")
        print(f"  Contact DOFs: {len(contact_params.get('contact_dofs', []))}")

        # Assemble frequency-domain matrix
        Z = self.assemble_frequency_matrix(omega)

        # Flatten external force
        F_ext_flat = F_ext.flatten()

        # Initialize with linear solution (no contact)
        print("\n  Step 1: Computing linear solution (initial guess)...")
        U_freq_flat = solve(Z, F_ext_flat)
        U_freq = U_freq_flat.reshape(self.n_dof, -1)

        # AFT iteration
        print(f"\n  Step 2: AFT Iteration (max {max_iter} iterations, tol={tol:.1e})...")

        history = {
            'residual': [],
            'max_force': [],
            'iteration': []
        }

        for iteration in range(max_iter):
            # Store previous solution
            U_freq_old = U_freq.copy()

            # 1. Frequency → Time
            u_time, t = self.frequency_to_time(U_freq, omega)

            # 2. Evaluate nonlinear contact forces in time domain
            f_nl_time = self.evaluate_contact_forces(u_time, contact_params)

            # 3. Time → Frequency
            F_nl_freq = self.time_to_frequency(f_nl_time, omega)
            F_nl_flat = F_nl_freq.flatten()

            # 4. Solve linear system with nonlinear forcing
            #    Z*U = F_ext - F_nl
            U_freq_new_flat = solve(Z, F_ext_flat - F_nl_flat)
            U_freq_new = U_freq_new_flat.reshape(self.n_dof, -1)

            # 5. Relaxation for stability
            U_freq = relaxation * U_freq_new + (1 - relaxation) * U_freq_old

            # 6. Check convergence
            residual = np.linalg.norm(U_freq - U_freq_old) / (np.linalg.norm(U_freq) + 1e-12)
            max_contact_force = np.max(np.abs(f_nl_time))

            history['residual'].append(residual)
            history['max_force'].append(max_contact_force)
            history['iteration'].append(iteration + 1)

            if (iteration + 1) % 10 == 0 or iteration < 5:
                print(f"    Iter {iteration+1:3d}: residual = {residual:.3e}, "
                      f"max contact force = {max_contact_force:.3e} N")

            if residual < tol:
                print(f"\n  ✓ Converged in {iteration+1} iterations!")
                break
        else:
            warnings.warn(f"AFT did not converge in {max_iter} iterations. "
                         f"Final residual: {residual:.3e}")

        # Final time-domain solution
        u_time, t = self.frequency_to_time(U_freq, omega)

        # Statistics
        f_nl_time_final = self.evaluate_contact_forces(u_time, contact_params)

        print(f"\n  Solution Statistics:")
        print(f"    Max displacement: {np.max(np.abs(u_time)):.3e} m")
        print(f"    Max contact force: {np.max(np.abs(f_nl_time_final)):.3e} N")
        print(f"    Contact active: {np.sum(f_nl_time_final > 1e-6)} / {f_nl_time_final.size} time-DOF pairs")

        return U_freq, u_time, history


def example_usage():
    """Demonstrate Harmonic Balance AFT solver."""
    print("Example: Harmonic Balance with Contact Nonlinearity")
    print("="*70)

    # Simple 2-DOF system with contact
    n_dof = 2

    # Mass and stiffness (arbitrary units)
    M = np.array([[1.0, 0.0],
                  [0.0, 0.5]])

    K = np.array([[1000.0, -500.0],
                  [-500.0, 1000.0]])

    # Initialize solver
    hb = HarmonicBalanceAFT(K, M, damping_ratio=0.02)
    hb.set_harmonics(n_harmonics=3)  # Use 3 harmonics

    # Excitation frequency and force
    omega = 10.0  # rad/s

    # External force (harmonic at fundamental frequency)
    F_ext = np.zeros((n_dof, 2*hb.n_harmonics + 1))
    F_ext[0, 1] = 100.0  # Cosine component at DOF 0

    # Contact parameters (unilateral contact at DOF 1)
    contact_params = {
        'contact_dofs': [1],
        'contact_stiffness': 5000.0,
        'gap_initial': 0.001,  # 1 mm initial gap
        'contact_type': 'penalty'
    }

    # Solve
    U_freq, u_time, history = hb.solve_harmonic_balance(
        omega, F_ext, contact_params,
        max_iter=50, tol=1e-6, relaxation=0.5
    )

    # Plot results
    try:
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(12, 8))

        T = 2 * np.pi / omega
        t = np.linspace(0, T, hb.n_time_points, endpoint=False)

        # Time-domain displacement
        axes[0, 0].plot(t*1000, u_time[0, :]*1000, 'b-', label='DOF 0')
        axes[0, 0].plot(t*1000, u_time[1, :]*1000, 'r-', label='DOF 1 (contact)')
        axes[0, 0].axhline(contact_params['gap_initial']*1000,
                          color='k', linestyle='--', label='Initial gap')
        axes[0, 0].set_xlabel('Time (ms)')
        axes[0, 0].set_ylabel('Displacement (mm)')
        axes[0, 0].set_title('Time-Domain Response')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)

        # Contact force
        f_nl_time = hb.evaluate_contact_forces(u_time, contact_params)
        axes[0, 1].plot(t*1000, f_nl_time[1, :], 'r-')
        axes[0, 1].set_xlabel('Time (ms)')
        axes[0, 1].set_ylabel('Contact Force (N)')
        axes[0, 1].set_title('Nonlinear Contact Force')
        axes[0, 1].grid(True, alpha=0.3)

        # Convergence history
        axes[1, 0].semilogy(history['iteration'], history['residual'], 'b-o')
        axes[1, 0].set_xlabel('Iteration')
        axes[1, 0].set_ylabel('Residual')
        axes[1, 0].set_title('Convergence History')
        axes[1, 0].grid(True, alpha=0.3)

        # Frequency content
        harmonics = range(hb.n_harmonics + 1)
        amplitudes_0 = [np.sqrt(U_freq[0, 0]**2)] + \
                      [np.sqrt(U_freq[0, 2*k-1]**2 + U_freq[0, 2*k]**2)
                       for k in range(1, hb.n_harmonics+1)]
        amplitudes_1 = [np.sqrt(U_freq[1, 0]**2)] + \
                      [np.sqrt(U_freq[1, 2*k-1]**2 + U_freq[1, 2*k]**2)
                       for k in range(1, hb.n_harmonics+1)]

        axes[1, 1].bar(np.array(harmonics)-0.15, amplitudes_0, 0.3,
                      label='DOF 0', alpha=0.7)
        axes[1, 1].bar(np.array(harmonics)+0.15, amplitudes_1, 0.3,
                      label='DOF 1 (contact)', alpha=0.7)
        axes[1, 1].set_xlabel('Harmonic Number')
        axes[1, 1].set_ylabel('Amplitude')
        axes[1, 1].set_title('Frequency Content')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('examples/output/harmonic_balance_example.png', dpi=150)
        print(f"\n✓ Example plot saved: examples/output/harmonic_balance_example.png")

    except ImportError:
        print("\n⚠️  Matplotlib not available for plotting")


if __name__ == "__main__":
    example_usage()
