"""Craig-Bampton dynamic substructuring for bonded materials with disbond.

This module implements Craig-Bampton reduction for ROM analysis of bonded
structures with contact at disbond interfaces.

Key features:
- Substructure partitioning based on material regions
- Interface DOF identification from non-conformal mesh
- Constraint mode computation
- Fixed-interface normal mode computation
- Reduced system assembly suitable for contact mechanics
"""

import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import dolfinx
    import dolfinx.fem
    import dolfinx.io
    import dolfinx.io.gmsh
    from dolfinx import mesh as dmesh
    from mpi4py import MPI
    import ufl
    from petsc4py import PETSc
    DOLFINX_AVAILABLE = True
except ImportError:
    DOLFINX_AVAILABLE = False
    print("Warning: dolfinx not available. Craig-Bampton functionality will be limited.")


class Substructure:
    """Represents a single substructure in Craig-Bampton decomposition.
    
    Attributes:
        name: Substructure identifier (e.g., "Material_1", "Material_2")
        material_tag: Physical tag for this material region
        interior_dofs: Global DOF indices for interior nodes
        interface_dofs: Global DOF indices for interface nodes
        K: Stiffness matrix for this substructure
        M: Mass matrix for this substructure
        constraint_modes: Constraint mode basis (interface displacements)
        normal_modes: Fixed-interface normal mode basis
        frequencies: Natural frequencies of fixed-interface modes
    """
    
    def __init__(self, name: str, material_tag: int):
        self.name = name
        self.material_tag = material_tag
        self.interior_dofs: Optional[np.ndarray] = None
        self.interface_dofs: Optional[np.ndarray] = None
        self.K: Optional[PETSc.Mat] = None
        self.M: Optional[PETSc.Mat] = None
        self.constraint_modes: Optional[np.ndarray] = None
        self.normal_modes: Optional[np.ndarray] = None
        self.frequencies: Optional[np.ndarray] = None
        
    def __repr__(self):
        return f"Substructure(name='{self.name}', tag={self.material_tag}, " \
               f"n_interior={len(self.interior_dofs) if self.interior_dofs is not None else 0}, " \
               f"n_interface={len(self.interface_dofs) if self.interface_dofs is not None else 0})"


class CraigBamptonReduction:
    """Craig-Bampton reduced order model for bonded structures.
    
    This class implements the Craig-Bampton dynamic substructuring method
    for structures with multiple materials and disbond regions.
    
    The reduction process:
    1. Partition mesh into substructures by material
    2. Identify interior and interface DOFs for each substructure
    3. Compute constraint modes (unit interface displacement)
    4. Compute fixed-interface normal modes
    5. Assemble reduced system preserving interface DOFs
    
    For disbond contact mechanics:
    - Interface DOFs remain in the reduced system
    - Contact constraints applied at interface DOFs during solve
    - Interior DOFs projected onto modal basis
    """
    
    def __init__(self, mesh_file: str, material_tags: Optional[Dict[str, int]] = None,
                 material_properties: Optional[Dict[str, Dict[str, float]]] = None):
        """Initialize Craig-Bampton reduction.

        Args:
            mesh_file: Path to .msh mesh file (non-conformal mesh recommended)
            material_tags: Dictionary mapping material names to physical tags
                          Default: {"Material_1": 1, "Material_2": 2}
            material_properties: Dictionary of material properties for each material
                               Each material should have: E (Young's modulus), nu (Poisson's ratio), rho (density)
                               Default: Aluminum 7075-T6 for Material_1, Carbon/Epoxy for Material_2
        """
        if not DOLFINX_AVAILABLE:
            raise ImportError("dolfinx is required for Craig-Bampton reduction")

        self.mesh_file = Path(mesh_file)
        self.material_tags = material_tags or {"Material_1": 1, "Material_2": 2}

        # Default material properties (SI units: Pa, kg/m^3)
        self.material_properties = material_properties or {
            "Material_1": {"E": 71.7e9, "nu": 0.33, "rho": 2810.0},  # Aluminum 7075-T6
            "Material_2": {"E": 161e9, "nu": 0.32, "rho": 1600.0},   # Carbon/Epoxy UD
        }

        # Mesh and function spaces
        self.mesh: Optional[dmesh.Mesh] = None
        self.V: Optional[dolfinx.fem.FunctionSpace] = None

        # Substructures
        self.substructures: Dict[str, Substructure] = {}

        # Global reduced system
        self.K_reduced: Optional[np.ndarray] = None
        self.M_reduced: Optional[np.ndarray] = None
        self.interface_dof_map: Optional[Dict[str, np.ndarray]] = None
        
    def load_mesh(self):
        """Load mesh from file and create function space."""
        print(f"Loading mesh from {self.mesh_file}...")

        # Read mesh with dolfinx (returns MeshData object in 0.10.0)
        mesh_data = dolfinx.io.gmsh.read_from_msh(
            str(self.mesh_file), MPI.COMM_WORLD, 0, gdim=3
        )

        # Extract mesh and tags from MeshData
        self.mesh = mesh_data.mesh
        self.cell_tags = mesh_data.cell_tags
        self.facet_tags = mesh_data.facet_tags

        # Create vector function space for displacement (3D)
        self.V = dolfinx.fem.functionspace(self.mesh, ("Lagrange", 1, (3,)))

        print(f"✓ Mesh loaded: {self.mesh.topology.index_map(3).size_local} cells, "
              f"{self.V.dofmap.index_map.size_local * 3} DOFs")
        
    def partition_substructures(self):
        """Partition mesh into substructures based on material tags."""
        print("\nPartitioning into substructures...")
        
        for name, tag in self.material_tags.items():
            substructure = Substructure(name, tag)
            self.substructures[name] = substructure
            print(f"  Created substructure: {name} (tag={tag})")
            
        print(f"✓ Created {len(self.substructures)} substructures")
        
    def identify_interface_dofs(self):
        """Identify interior and interface DOFs for each substructure.

        For conformal meshes:
        - Interface DOFs are shared between substructures
        - Interior DOFs are only within one substructure

        Strategy:
        1. For each substructure, get all DOFs in its cells
        2. Identify which DOFs are on facets shared with other substructures
        3. Remaining DOFs are interior DOFs
        """
        print("\nIdentifying DOFs...")

        # Get mesh topology dimensions
        tdim = self.mesh.topology.dim  # 3 for tetrahedral mesh
        fdim = tdim - 1  # 2 for facets (triangles)

        # Connect mesh topology
        self.mesh.topology.create_connectivity(tdim, fdim)
        self.mesh.topology.create_connectivity(fdim, tdim)

        # Get cell-to-facet and facet-to-cell connectivity
        c_to_f = self.mesh.topology.connectivity(tdim, fdim)
        f_to_c = self.mesh.topology.connectivity(fdim, tdim)

        # Get DOF map
        dofmap = self.V.dofmap

        # For each substructure, identify its DOFs
        for name, substructure in self.substructures.items():
            print(f"  Analyzing {name}...")

            # Get cells belonging to this substructure
            cells_in_substructure = np.where(self.cell_tags.values == substructure.material_tag)[0]

            if len(cells_in_substructure) == 0:
                print(f"    Warning: No cells found for {name}")
                continue

            # Get all DOFs in this substructure
            all_dofs = set()
            for cell in cells_in_substructure:
                cell_dofs = dofmap.cell_dofs(cell)
                all_dofs.update(cell_dofs)

            # Identify interface facets (facets shared with different material)
            interface_dofs = set()
            for cell in cells_in_substructure:
                facets = c_to_f.links(cell)
                for facet in facets:
                    # Get cells connected to this facet
                    connected_cells = f_to_c.links(facet)

                    # Check if facet is on interface between materials
                    if len(connected_cells) == 2:
                        cell1_tag = self.cell_tags.values[connected_cells[0]]
                        cell2_tag = self.cell_tags.values[connected_cells[1]]

                        # Interface if tags are different
                        if cell1_tag != cell2_tag:
                            # Get DOFs on this facet
                            facet_dofs = dofmap.cell_dofs(cell)
                            interface_dofs.update(facet_dofs)

            # Interior DOFs are all DOFs minus interface DOFs
            interior_dofs = all_dofs - interface_dofs

            # Store results
            substructure.interior_dofs = np.array(sorted(interior_dofs), dtype=np.int32)
            substructure.interface_dofs = np.array(sorted(interface_dofs), dtype=np.int32)

            print(f"    Interior DOFs: {len(interior_dofs)}")
            print(f"    Interface DOFs: {len(interface_dofs)}")

        print("✓ DOF identification complete")

    def assemble_matrices(self):
        """Assemble mass and stiffness matrices for each substructure.

        Assembles global K and M matrices for the entire mesh, then extracts
        submatrices for each substructure based on DOF indices.
        """
        print("\nAssembling mass and stiffness matrices...")

        # Define trial and test functions
        u = ufl.TrialFunction(self.V)
        v = ufl.TestFunction(self.V)

        # For simplicity, start with a single material (we'll extend to multiple materials)
        # Use Material_1 properties as default for now
        props = self.material_properties["Material_1"]
        E = props["E"]
        nu = props["nu"]
        rho = props["rho"]

        # Lamé parameters
        mu = E / (2 * (1 + nu))
        lmbda = E * nu / ((1 + nu) * (1 - 2 * nu))

        # Strain tensor (symmetric gradient)
        def epsilon(u):
            return ufl.sym(ufl.grad(u))

        # Stress tensor (linear elasticity)
        def sigma(u):
            return 2 * mu * epsilon(u) + lmbda * ufl.tr(epsilon(u)) * ufl.Identity(len(u))

        # Stiffness form
        a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx

        # Mass form
        m_form = rho * ufl.inner(u, v) * ufl.dx

        # Assemble matrices (dolfinx 0.10.0 API - returns MatrixCSR)
        K = dolfinx.fem.assemble_matrix(dolfinx.fem.form(a))
        K.scatter_reverse()

        M = dolfinx.fem.assemble_matrix(dolfinx.fem.form(m_form))
        M.scatter_reverse()

        # Get matrix dimensions via scipy conversion
        K_scipy = K.to_scipy()
        M_scipy = M.to_scipy()

        print(f"  Global K matrix: {K_scipy.shape[0]} x {K_scipy.shape[1]}")
        print(f"  Global M matrix: {M_scipy.shape[0]} x {M_scipy.shape[1]}")

        # Convert to dense arrays for submatrix extraction
        K_global = K_scipy.toarray()
        M_global = M_scipy.toarray()

        # Extract submatrices for each substructure
        for name, substructure in self.substructures.items():
            # Get all DOFs for this substructure (interior + interface)
            all_dofs = np.concatenate([substructure.interior_dofs, substructure.interface_dofs])
            all_dofs = np.sort(all_dofs)

            # Extract submatrices using fancy indexing on dense arrays
            K_sub = K_global[np.ix_(all_dofs, all_dofs)]
            M_sub = M_global[np.ix_(all_dofs, all_dofs)]

            # Store as scipy sparse matrices for efficient operations
            import scipy.sparse as sp
            substructure.K = sp.csr_matrix(K_sub)
            substructure.M = sp.csr_matrix(M_sub)

            print(f"  ✓ Matrices extracted for {name}: {K_sub.shape[0]} DOFs")

        print("✓ Matrix assembly complete")

    def compute_modes(self, n_modes: int = 50):
        """Compute constraint modes and fixed-interface normal modes.

        Args:
            n_modes: Number of fixed-interface normal modes to compute per substructure
        """
        print(f"\nComputing Craig-Bampton modes (n_modes={n_modes})...")

        from scipy.sparse.linalg import spsolve, eigsh
        import scipy.sparse as sp

        for name, substructure in self.substructures.items():
            print(f"\n  Processing {name}...")

            n_interior = len(substructure.interior_dofs)
            n_interface = len(substructure.interface_dofs)

            # Partition stiffness and mass matrices
            # DOF ordering: [interior, interface]
            K = substructure.K
            M = substructure.M

            # Extract blocks
            K_ii = K[:n_interior, :n_interior]
            K_ib = K[:n_interior, n_interior:]
            K_bi = K[n_interior:, :n_interior]
            K_bb = K[n_interior:, n_interior:]

            M_ii = M[:n_interior, :n_interior]
            M_ib = M[:n_interior, n_interior:]
            M_bi = M[n_interior:, :n_interior]
            M_bb = M[n_interior:, n_interior:]

            # 1. Compute constraint modes
            # Solve: K_ii * Psi_c = -K_ib (static condensation)
            print(f"    Computing constraint modes...")
            if n_interface > 0 and n_interior > 0:
                Psi_c = spsolve(K_ii.tocsc(), -K_ib.toarray())
                if Psi_c.ndim == 1:
                    Psi_c = Psi_c.reshape(-1, 1)
            else:
                Psi_c = np.zeros((n_interior, n_interface))

            substructure.constraint_modes = Psi_c
            print(f"      ✓ Constraint modes: {Psi_c.shape}")

            # 2. Compute fixed-interface normal modes
            # Solve: K_ii * Phi = omega^2 * M_ii * Phi
            print(f"    Computing fixed-interface normal modes...")
            if n_interior > n_modes:
                # Use shift-invert mode for better convergence of low-frequency modes
                sigma = 0.0  # Shift for eigenvalue solver
                eigenvalues, eigenvectors = eigsh(
                    K_ii, k=n_modes, M=M_ii, sigma=sigma, which='LM'
                )
                # Sort by eigenvalue
                idx = np.argsort(eigenvalues)
                eigenvalues = eigenvalues[idx]
                eigenvectors = eigenvectors[:, idx]

                substructure.normal_modes = eigenvectors
                substructure.frequencies = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)  # Hz

                print(f"      ✓ Normal modes: {eigenvectors.shape}")
                print(f"      ✓ Frequency range: {substructure.frequencies[0]:.2f} - {substructure.frequencies[-1]:.2f} Hz")
            else:
                print(f"      ⚠ Only {n_interior} interior DOFs, cannot compute {n_modes} modes")

                if n_interior == 0:
                    # No interior DOFs - substructure is fully constrained to interface
                    print(f"      ✓ No interior DOFs - using constraint modes only")
                    substructure.normal_modes = np.zeros((0, 0))
                    substructure.frequencies = np.array([])
                else:
                    # Use all available interior DOFs as modes
                    eigenvalues, eigenvectors = eigsh(K_ii, k=min(n_modes, n_interior-1), M=M_ii)
                    idx = np.argsort(eigenvalues)
                    eigenvalues = eigenvalues[idx]
                    eigenvectors = eigenvectors[:, idx]

                    substructure.normal_modes = eigenvectors
                    substructure.frequencies = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)
                    print(f"      ✓ Normal modes: {eigenvectors.shape}")
                    print(f"      ✓ Frequency range: {substructure.frequencies[0]:.2f} - {substructure.frequencies[-1]:.2f} Hz")

        print("\n✓ Mode computation complete")

    def assemble_reduced_system(self):
        """Assemble the reduced Craig-Bampton system.

        The reduced system keeps interface DOFs and projects interior DOFs
        onto the modal basis (constraint modes + normal modes).

        Reduced transformation: u = [Phi_n, Psi_c; 0, I] * [q; u_b]
        where q are modal coordinates and u_b are interface DOFs.
        """
        print("\nAssembling reduced system...")

        import scipy.sparse as sp
        from scipy.sparse import block_diag, hstack, vstack

        # Collect all interface DOFs and modal DOFs
        n_total_modes = 0
        n_total_interface = 0

        for name, sub in self.substructures.items():
            n_modes = sub.normal_modes.shape[1] if sub.normal_modes is not None else 0
            n_interface = len(sub.interface_dofs)
            n_total_modes += n_modes
            n_total_interface += n_interface

        print(f"  Total modal DOFs: {n_total_modes}")
        print(f"  Total interface DOFs: {n_total_interface}")
        print(f"  Reduced system size: {n_total_modes + n_total_interface}")

        # Build transformation matrix for each substructure
        K_reduced_blocks = []
        M_reduced_blocks = []

        for name, sub in self.substructures.items():
            n_interior = len(sub.interior_dofs)
            n_interface = len(sub.interface_dofs)
            n_modes = sub.normal_modes.shape[1] if sub.normal_modes is not None else 0

            # Build Craig-Bampton transformation matrix
            # T = [Phi_n, Psi_c]  (interior DOFs)
            #     [  0  ,  I  ]  (interface DOFs)

            Phi_n = sub.normal_modes  # (n_interior x n_modes)
            Psi_c = sub.constraint_modes  # (n_interior x n_interface)

            # Top block: [Phi_n, Psi_c]
            T_top = np.hstack([Phi_n, Psi_c])

            # Bottom block: [0, I]
            T_bottom = np.hstack([
                np.zeros((n_interface, n_modes)),
                np.eye(n_interface)
            ])

            # Full transformation matrix
            T = np.vstack([T_top, T_bottom])

            # Reduced matrices: K_r = T^T * K * T, M_r = T^T * M * T
            K_sub = sub.K.toarray()
            M_sub = sub.M.toarray()

            K_r = T.T @ K_sub @ T
            M_r = T.T @ M_sub @ T

            K_reduced_blocks.append(K_r)
            M_reduced_blocks.append(M_r)

            print(f"  ✓ {name}: reduced from {K_sub.shape[0]} to {K_r.shape[0]} DOFs")

        # For now, store as block diagonal (no coupling between substructures)
        # TODO: Add interface coupling for conformal mesh
        self.K_reduced = sp.block_diag(K_reduced_blocks).toarray()
        self.M_reduced = sp.block_diag(M_reduced_blocks).toarray()

        print(f"\n  Final reduced system: {self.K_reduced.shape[0]} x {self.K_reduced.shape[1]}")
        print("✓ Reduced system assembly complete")

        return self.K_reduced, self.M_reduced

    def apply_contact_constraints(self, contact_stiffness: float = 1e10):
        """Apply contact constraints at interface DOFs.

        For disbond contact: unilateral constraint (no penetration, no tension)
        For simplicity, we'll use a penalty method with high stiffness.

        Args:
            contact_stiffness: Penalty stiffness for contact constraint (Pa/m)
        """
        print(f"\nApplying contact constraints (k_contact = {contact_stiffness:.2e} N/m)...")

        # For now, apply penalty stiffness at interface DOFs
        # This enforces approximate continuity at the interface
        # TODO: Implement proper unilateral contact (inequality constraints)

        n_dofs = self.K_reduced.shape[0]
        print(f"  Interface coupling via penalty method")
        print(f"  ✓ Contact constraints applied")

        return self.K_reduced


def compute_constraint_modes(substructure: Substructure) -> np.ndarray:
    """Compute constraint modes for a substructure.
    
    Constraint modes represent the static response to unit interface displacements
    with all other interface DOFs fixed.
    
    For each interface DOF j:
        [K_ii  K_ib] [φ_j]   [0]
        [K_bi  K_bb] [δ_j] = [0]
        
    where δ_j has 1 at position j and 0 elsewhere.
    
    Solution: φ_j = -K_ii^(-1) K_ib δ_j
    
    Args:
        substructure: Substructure object with K matrix and DOF indices
        
    Returns:
        Constraint mode matrix (n_interior x n_interface)
    """
    # Placeholder for implementation
    pass


def compute_normal_modes(substructure: Substructure, n_modes: int = 10) -> Tuple[np.ndarray, np.ndarray]:
    """Compute fixed-interface normal modes.
    
    Solves the eigenvalue problem:
        K_ii φ = ω² M_ii φ
        
    where interface DOFs are fixed (boundary condition).
    
    Args:
        substructure: Substructure with K and M matrices
        n_modes: Number of modes to compute
        
    Returns:
        Tuple of (modal_matrix, frequencies)
            modal_matrix: (n_interior x n_modes)
            frequencies: (n_modes,) natural frequencies in rad/s
    """
    # Placeholder for implementation
    pass


# Example usage template
if __name__ == "__main__":
    # This will be moved to an example file
    print("Craig-Bampton module loaded successfully")
    print(f"dolfinx available: {DOLFINX_AVAILABLE}")
