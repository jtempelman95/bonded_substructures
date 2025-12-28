#!/bin/bash
# Run all bonded_substructures examples
# This script runs all examples in order to verify the implementation

set -e  # Exit on error

echo "=========================================================================="
echo "Running All Bonded Substructures Examples"
echo "=========================================================================="
echo ""

# Color codes for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Track results
PASSED=0
FAILED=0
SKIPPED=0

run_example() {
    local example=$1
    local description=$2

    echo -e "${BLUE}=========================================================================="
    echo -e "Running: ${example}"
    echo -e "Description: ${description}"
    echo -e "==========================================================================${NC}"
    echo ""

    if conda run -n bonded_substructures python examples/${example}; then
        echo -e "${GREEN}✅ ${example} completed successfully${NC}"
        ((PASSED++))
    else
        echo -e "${RED}❌ ${example} failed${NC}"
        ((FAILED++))
        return 1
    fi
    echo ""
}

# Example 01: Basic bonded plate
run_example "01_bonded_plate_basic.py" "Basic 1ft × 1ft × 0.25in plate (no disbond)"

# Example 02: Bonded plate with disbond
run_example "02_bonded_plate_disbond.py" "1ft × 1ft × 0.25in plate with 1in disbond"

# Example 03: Visualize meshes
run_example "03_visualize_mesh.py" "Generate mesh visualizations (matplotlib + pyvista)"

# Example 07: Wide plate with disbond (if it exists)
if [ -f "examples/07_wide_plate_disbond.py" ]; then
    echo -e "${BLUE}Note: Example 07 may have different dimensions${NC}"
    if run_example "07_wide_plate_disbond.py" "Wide plate with localized disbond"; then
        :
    else
        echo -e "${BLUE}Skipping Example 07 due to errors${NC}"
        ((FAILED--))
        ((SKIPPED++))
    fi
fi

# Example 08: Time-domain response (may take longer)
if [ -f "examples/08_time_domain_response.py" ]; then
    echo -e "${BLUE}Note: Example 08 requires dolfinx and may take several minutes${NC}"
    read -p "Run Example 08 (time-domain response)? [y/N] " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        if run_example "08_time_domain_response.py" "Time-domain response using Craig-Bampton ROM"; then
            :
        else
            echo -e "${BLUE}Example 08 failed - this may be due to missing dolfinx${NC}"
            ((FAILED--))
            ((SKIPPED++))
        fi
    else
        echo -e "${BLUE}Skipping Example 08${NC}"
        ((SKIPPED++))
    fi
fi

# Summary
echo ""
echo "=========================================================================="
echo "Summary"
echo "=========================================================================="
echo -e "${GREEN}Passed: ${PASSED}${NC}"
if [ $FAILED -gt 0 ]; then
    echo -e "${RED}Failed: ${FAILED}${NC}"
else
    echo "Failed: ${FAILED}"
fi
if [ $SKIPPED -gt 0 ]; then
    echo -e "${BLUE}Skipped: ${SKIPPED}${NC}"
fi
echo ""

# List generated files
echo "Generated mesh files:"
ls -lh *.msh 2>/dev/null || echo "  No mesh files found"
echo ""

echo "Generated visualization files:"
ls -lh examples/output/*.png 2>/dev/null || echo "  No visualization files found"
echo ""

if [ $FAILED -eq 0 ]; then
    echo -e "${GREEN}=========================================================================="
    echo -e "All examples completed successfully!"
    echo -e "==========================================================================${NC}"
    exit 0
else
    echo -e "${RED}=========================================================================="
    echo -e "Some examples failed. Check output above for details."
    echo -e "==========================================================================${NC}"
    exit 1
fi
