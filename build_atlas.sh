#!/bin/bash
#
# Build Atlas and dependencies from source for FESOM2 tracer dwarf
# using Atlas's supported tools/install.sh installer.
#
# Usage:
#   ./build_atlas.sh [gnu|intel|nvidia] [--clean] [--verbose]
#   ./build_atlas.sh --compiler gnu|intel|nvidia [--clean] [--verbose]
#

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
COMPILER="gnu"
VERBOSE=false
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
            VERBOSE=true
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

# Installation prefix (all dependencies + Atlas go here)
INSTALL_PREFIX="${SCRIPT_DIR}/atlas_install_${COMPILER}"
WORK_DIR="${SCRIPT_DIR}/atlas_builds_${COMPILER}"
ATLAS_SOURCE="${SCRIPT_DIR}/atlas_deps/atlas"
ATLAS_INSTALLER="${ATLAS_SOURCE}/tools/install.sh"
ATLAS_GIT_REPOSITORY="${ATLAS_GIT_REPOSITORY:-https://github.com/ecmwf/atlas.git}"
ATLAS_GIT_VERSION="${ATLAS_GIT_VERSION:-develop}"

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

if ! command -v "$FC" >/dev/null 2>&1; then
    echo "Error: Fortran compiler not found: $FC"
    exit 1
fi
if ! command -v "$CC" >/dev/null 2>&1; then
    echo "Error: C compiler not found: $CC"
    exit 1
fi
if ! command -v "$CXX" >/dev/null 2>&1; then
    echo "Error: C++ compiler not found: $CXX"
    exit 1
fi

if command -v nproc >/dev/null 2>&1; then
    JOBS=$(nproc)
else
    JOBS=$(sysctl -n hw.ncpu)
fi

echo "========================================="
echo "Atlas Dependency Build"
echo "========================================="
echo "Compiler:  $COMPILER"
echo "FC:        $FC"
echo "CC:        $CC"
echo "CXX:       $CXX"
echo "Install:   $INSTALL_PREFIX"
echo "Work dir:  $WORK_DIR"
echo "========================================="
echo ""

# Clean if requested
if $DO_CLEAN; then
    echo "Cleaning old builds..."
    rm -rf "$INSTALL_PREFIX" "$WORK_DIR"
fi

if [ ! -d "${ATLAS_SOURCE}/.git" ]; then
    if [ -e "$ATLAS_SOURCE" ]; then
        echo "Error: Atlas source exists but is not a Git checkout: $ATLAS_SOURCE"
        exit 1
    fi
    echo "Cloning Atlas $ATLAS_GIT_VERSION..."
    mkdir -p "$(dirname "$ATLAS_SOURCE")"
    git clone --depth 1 --branch "$ATLAS_GIT_VERSION" \
        "$ATLAS_GIT_REPOSITORY" "$ATLAS_SOURCE"
fi

if [ ! -f "$ATLAS_INSTALLER" ]; then
    echo "Error: Atlas installer not found: $ATLAS_INSTALLER"
    exit 1
fi

mkdir -p "$WORK_DIR"
FCKIT_CACHE="${WORK_DIR}/builds/fckit/CMakeCache.txt"
if [ -f "$FCKIT_CACHE" ]; then
    CACHED_ECKIT_DIR=$(sed -n 's/^eckit_DIR:PATH=//p' "$FCKIT_CACHE")
    case "$CACHED_ECKIT_DIR" in
        ""|"${INSTALL_PREFIX}"/*) ;;
        *)
            echo "Removing stale fckit cache referencing: $CACHED_ECKIT_DIR"
            rm -rf "${WORK_DIR}/builds/fckit"
            ;;
    esac
fi

CMAKE_OPTIONS="-DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_Fortran_COMPILER=$FC -DENABLE_MPI=ON -DENABLE_OMP=OFF"
if $VERBOSE; then
    CMAKE_OPTIONS="$CMAKE_OPTIONS -DCMAKE_VERBOSE_MAKEFILE=ON"
fi

# Do not let package hints from another Atlas installation override the
# dependencies that tools/install.sh installs into INSTALL_PREFIX.
unset ecbuild_DIR eckit_DIR fckit_DIR atlas_DIR
unset ecbuild_ROOT eckit_ROOT fckit_ROOT atlas_ROOT
unset ECBUILD_ROOT ECKIT_ROOT FCKIT_ROOT ATLAS_ROOT
export CC CXX FC
bash "$ATLAS_INSTALLER" \
    --with-deps \
    --enable-fortran \
    --enable-lz4 \
    --with-atlas-fesom \
    --prefix "$INSTALL_PREFIX" \
    --work-dir "$WORK_DIR" \
    --build-type Release \
    --parallel "$JOBS" \
    --cmake "$CMAKE_OPTIONS"

echo ""
echo "========================================="
echo "Build complete!"
echo "========================================="
echo "Install prefix: $INSTALL_PREFIX"
echo "Export this for use in configure.sh:"
echo "  export atlas_DIR=$INSTALL_PREFIX/lib/cmake/atlas"
echo ""
