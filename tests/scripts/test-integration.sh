#!/bin/bash
# Integration Test Runner for PEPPRO
# Uses bulker exec to run pytest with bioinformatics tools available
# via containerized commands.
#
# Usage:
#   ./tests/scripts/test-integration.sh                  # Run all integration tests
#   ./tests/scripts/test-integration.sh -k "test_se"     # Run specific tests
#   ./tests/scripts/test-integration.sh --keep-test-outputs  # Preserve outputs

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
TESTS_DIR="$SCRIPT_DIR/.."

BULKER_CRATE="${PEPPRO_TEST_BULKER_CRATE:-databio/peppro:1.1.0}"

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${GREEN}=== PEPPRO Integration Tests ===${NC}"
echo -e "\n${GREEN}Running integration tests via bulker exec ${BULKER_CRATE}...${NC}"
echo ""

cd "$PROJECT_ROOT"

set +e
RUN_INTEGRATION_TESTS=true bulker exec "${BULKER_CRATE}" -- \
    python3 -m pytest "$TESTS_DIR/test_integration.py" -v "$@"
PYTEST_EXIT=$?
set -e

if [ $PYTEST_EXIT -eq 0 ]; then
    echo -e "\n${GREEN}Integration tests completed successfully!${NC}"
else
    echo -e "\n${RED}Integration tests failed (exit code: ${PYTEST_EXIT})${NC}"
fi
exit $PYTEST_EXIT
