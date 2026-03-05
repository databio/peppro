#!/bin/bash
# Integration Test Runner for PEPPRO
# Activates the bulker crate, then runs pytest integration tests with
# bioinformatics tools available via containerized commands.
#
# Usage:
#   ./tests/scripts/test-integration.sh                  # Run all integration tests
#   ./tests/scripts/test-integration.sh -k "test_se"     # Run specific tests
#   ./tests/scripts/test-integration.sh --keep-test-outputs  # Preserve outputs

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
TESTS_DIR="$SCRIPT_DIR/.."

BULKER_CRATE="${PEPPRO_TEST_BULKER_CRATE:-databio/peppro:1.0.12}"

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${GREEN}=== PEPPRO Integration Tests ===${NC}"

# Load crate from local manifest if not already loaded
MANIFEST="$TESTS_DIR/bulker_manifest.yaml"
if ! bulker list 2>/dev/null | grep -q "${BULKER_CRATE}"; then
    if [ -f "$MANIFEST" ]; then
        echo -e "${YELLOW}Crate ${BULKER_CRATE} not loaded. Loading from local manifest...${NC}"
        bulker load "${BULKER_CRATE}" -m "$MANIFEST" -r
    else
        echo -e "${RED}ERROR: Crate ${BULKER_CRATE} not loaded and no manifest at ${MANIFEST}${NC}"
        exit 1
    fi
fi

# Activate bulker crate by prepending its path
CRATE_PATH=$(bulker list 2>/dev/null | grep "${BULKER_CRATE}" | sed 's/.* -- //')
if [ -z "$CRATE_PATH" ]; then
    echo -e "${RED}ERROR: Could not find crate path for ${BULKER_CRATE}${NC}"
    exit 1
fi

export PATH="${CRATE_PATH}:${PATH}"
export RUN_INTEGRATION_TESTS=true
export PEPPRO_TEST_BULKER_CRATE="$BULKER_CRATE"

echo -e "\n${GREEN}Running integration tests...${NC}"
echo "  Crate PATH: ${CRATE_PATH}"
echo ""

cd "$PROJECT_ROOT"

set +e
python3 -m pytest "$TESTS_DIR/test_integration.py" -v "$@"
PYTEST_EXIT=$?
set -e

if [ $PYTEST_EXIT -eq 0 ]; then
    echo -e "\n${GREEN}Integration tests completed successfully!${NC}"
else
    echo -e "\n${RED}Integration tests failed (exit code: ${PYTEST_EXIT})${NC}"
fi
exit $PYTEST_EXIT
