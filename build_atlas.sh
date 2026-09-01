#!/bin/bash
#
# Compatibility wrapper: Atlas installation is owned by CMakeLists.txt.
#
# Usage:
#   ./build_atlas.sh [gnu|intel|nvidia] [--clean] [--verbose]
#   ./build_atlas.sh --compiler gnu|intel|nvidia [--clean] [--verbose]
#

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
COMPILER="gnu"
DO_CLEAN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        gnu|intel|nvidia)
            COMPILER="$1"
            shift
            ;;
        --compiler)
            COMPILER="$2"
            shift 2
            ;;
        --clean)
            DO_CLEAN=true
            shift
            ;;
        --verbose)
            shift
            ;;
        *)
            echo "Error: unknown option '$1'"
            echo "Usage: ./build_atlas.sh [gnu|intel|nvidia] [--clean] [--verbose]"
            echo "   or: ./build_atlas.sh --compiler gnu|intel|nvidia [--clean] [--verbose]"
            exit 1
            ;;
    esac
done

case "$COMPILER" in
    gnu)
        FC="${FC:-$(command -v gfortran 2>/dev/null || true)}"
        CC="${CC:-$(command -v gcc 2>/dev/null || true)}"
        CXX="${CXX:-$(command -v g++ 2>/dev/null || true)}"
        ;;
    intel)
        FC="${FC:-$(command -v ifx 2>/dev/null || true)}"
        CC="${CC:-$(command -v icx 2>/dev/null || true)}"
        CXX="${CXX:-$(command -v icpx 2>/dev/null || true)}"
        ;;
    nvidia)
        FC="${FC:-$(command -v nvfortran 2>/dev/null || true)}"
        CC="${CC:-$(command -v nvc 2>/dev/null || true)}"
        CXX="${CXX:-$(command -v nvc++ 2>/dev/null || true)}"
        ;;
esac


CMAKE_ARGS=(
    -S "$SCRIPT_DIR"
    -B "$SCRIPT_DIR/build_${COMPILER}_dp_atlas"
    -DENABLE_ATLAS=ON
    -DCMAKE_Fortran_COMPILER="$FC"
    -DCMAKE_C_COMPILER="$CC"
)
if $DO_CLEAN; then
    CMAKE_ARGS+=(-DATLAS_CLEAN_INSTALL=ON)
fi
exec cmake "${CMAKE_ARGS[@]}"
