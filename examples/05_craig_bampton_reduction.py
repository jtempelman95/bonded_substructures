"""Example 05: Craig-Bampton reduction for disbond ROM.

This example demonstrates the Craig-Bampton dynamic substructuring method
applied to a bonded structure with disbond.

Workflow:
1. Load non-conformal mesh (from Example 04)
2. Partition into substructures (Material 1 and Material 2)
3. Identify interior and interface DOFs
4. Compute constraint modes for each substructure
5. Compute fixed-interface normal modes
6. Assemble reduced system

The reduced system preserves interface DOFs for contact mechanics while
reducing interior DOFs using modal basis.
"""

from pathlib import Path

try:
    from bonded_substructures.rom.craig_bampton import CraigBamptonReduction
    CB_AVAILABLE = True
except ImportError as e:
    CB_AVAILABLE = False
    print(f"Craig-Bampton functionality not available: {e}")
    print("Install dolfinx to use this example")


def main():
    """Demonstrate Craig-Bampton reduction workflow."""
    
    print("=" * 70)
    print("Craig-Bampton Reduction for Disbond ROM")
    print("=" * 70)
    
    if not CB_AVAILABLE:
        print("\n⚠️  Craig-Bampton reduction requires dolfinx")
        print("Install with: conda install -c conda-forge fenics-dolfinx")
        return
    
    # Use conformal mesh for now (non-conformal mesh via meshio has compatibility issues with dolfinx)
    # TODO: Implement non-conformal mesh using gmsh API directly
    mesh_file = "rectangle_disbond.msh"

    if not Path(mesh_file).exists():
        print(f"\n⚠️  Mesh file not found: {mesh_file}")
        print("Run Example 02 first to generate disbond mesh")
        return

    print(f"\nInput mesh: {mesh_file}")
    print("(Conformal mesh - will implement non-conformal approach later)")
    
    # Initialize Craig-Bampton reduction
    print("\n" + "-" * 70)
    print("Step 1: Initialize and Load Mesh")
    print("-" * 70)
    
    cb = CraigBamptonReduction(
        mesh_file=mesh_file,
        material_tags={"Material_1": 1, "Material_2": 2}
    )
    
    try:
        cb.load_mesh()
    except Exception as e:
        print(f"Error loading mesh: {e}")
        return
    
    # Partition into substructures
    print("\n" + "-" * 70)
    print("Step 2: Partition into Substructures")
    print("-" * 70)
    
    cb.partition_substructures()
    
    # Identify DOFs
    print("\n" + "-" * 70)
    print("Step 3: Identify Interior and Interface DOFs")
    print("-" * 70)

    cb.identify_interface_dofs()

    # Print substructure info
    print("\nSubstructure Summary:")
    for name, sub in cb.substructures.items():
        print(f"  {sub}")

    # Assemble matrices
    print("\n" + "-" * 70)
    print("Step 4: Assemble Mass and Stiffness Matrices")
    print("-" * 70)

    cb.assemble_matrices()

    # Compute modes
    print("\n" + "-" * 70)
    print("Step 5: Compute Craig-Bampton Modes")
    print("-" * 70)

    n_modes = 50
    cb.compute_modes(n_modes=n_modes)

    # Assemble reduced system
    print("\n" + "-" * 70)
    print("Step 6: Assemble Reduced System")
    print("-" * 70)

    K_r, M_r = cb.assemble_reduced_system()

    # Apply contact constraints
    print("\n" + "-" * 70)
    print("Step 7: Apply Contact Constraints")
    print("-" * 70)

    cb.apply_contact_constraints(contact_stiffness=1e10)

    # Summary
    print("\n" + "-" * 70)
    print("Summary")
    print("-" * 70)
    print(f"Original system: 690 DOFs")
    print(f"Reduced system: {K_r.shape[0]} DOFs")
    print(f"Reduction factor: {690/K_r.shape[0]:.1f}x")
    print(f"\nNext step: Run Example 06 for harmonic response simulation")
    
    print("\n" + "=" * 70)
    print("Craig-Bampton Reduction Initialized Successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
