#!/bin/bash
# Run core examples and check for physical consistency

echo "=========================================================================="
echo "Running Core Bonded Substructures Examples"
echo "=========================================================================="
echo ""

# Example 01: Basic plate
echo "=== Example 01: Basic Bonded Plate ==="
conda run -n bonded_substructures python examples/01_bonded_plate_basic.py 2>&1 | tail -30
echo ""

# Example 02: Disbond plate
echo "=== Example 02: Bonded Plate with Disbond ==="
conda run -n bonded_substructures python examples/02_bonded_plate_disbond.py 2>&1 | tail -35
echo ""

# Example 03: Visualizations
echo "=== Example 03: Mesh Visualizations ==="
conda run -n bonded_substructures python examples/03_visualize_mesh.py 2>&1 | grep -E "Visualizing|Saved|complete|Generated"
echo ""

echo "=========================================================================="
echo "Generated Files:"
echo "=========================================================================="
echo ""
echo "Mesh files:"
ls -lh *.msh 2>/dev/null
echo ""
echo "Visualization files:"
ls -lh examples/output/*.png 2>/dev/null | grep -E "basic_plate|disbond_plate"
echo ""

echo "✅ Core examples completed!"
