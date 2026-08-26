#!/bin/bash
#
# Build Atlas and dependencies from source for FESOM2 tracer dwarf
# Uses system compiler and MPICH; installs to local directory
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

# Source directories
DEPS_DIR="${SCRIPT_DIR}/atlas_deps"
ECBUILD_SOURCE="${DEPS_DIR}/ecbuild"
ECKIT_SOURCE="${DEPS_DIR}/eckit"
FCKIT_SOURCE="${DEPS_DIR}/fckit"
ATLAS_SOURCE="${DEPS_DIR}/atlas"
ATLAS_GIT_REPOSITORY="https://github.com/ecmwf/atlas.git"
ATLAS_GIT_VERSION="fix/fortran-unstructured-grid-by-id"

# Build directories
BUILD_BASE="${SCRIPT_DIR}/atlas_builds_${COMPILER}"
ECBUILD_BUILD="${BUILD_BASE}/ecbuild"
ECKIT_BUILD="${BUILD_BASE}/eckit"
FCKIT_BUILD="${BUILD_BASE}/fckit"
ATLAS_BUILD="${BUILD_BASE}/atlas"

# Compiler setup
case "$COMPILER" in
    gnu)
        FC=$(which gfortran)
        CC=$(which gcc)
        CXX=$(which g++)
        ;;
    intel)
        source /opt/intel/oneapi/setvars.sh --force > /dev/null 2>&1 || true
        FC=$(which ifx 2>/dev/null || echo "/opt/intel/oneapi/compiler/latest/bin/ifx")
        CC=$(which icx 2>/dev/null || echo "/opt/intel/oneapi/compiler/latest/bin/icx")
        CXX=$(which icpx 2>/dev/null || echo "/opt/intel/oneapi/compiler/latest/bin/icpx")
        ;;
    nvidia)
        FC="/opt/nvidia/hpc_sdk/Linux_x86_64/2025/compilers/bin/nvfortran"
        CC="/opt/nvidia/hpc_sdk/Linux_x86_64/2025/compilers/bin/nvc"
        CXX="/opt/nvidia/hpc_sdk/Linux_x86_64/2025/compilers/bin/nvc++"
        ;;
    *)
        echo "Error: unknown compiler '$COMPILER'"
        exit 1
        ;;
esac

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
echo "Atlas:     $ATLAS_GIT_VERSION"
echo "========================================="
echo ""

# Clean if requested
if $DO_CLEAN; then
    echo "Cleaning old builds..."
    rm -rf "$INSTALL_PREFIX" "$BUILD_BASE"
fi

mkdir -p "$DEPS_DIR" "$BUILD_BASE"

# ========================================
# 1. Clone dependencies if not present
# ========================================
if [ ! -d "$ECBUILD_SOURCE" ]; then
    echo "Cloning ecbuild..."
    git clone --depth 1 --branch master https://github.com/ecmwf/ecbuild.git "$ECBUILD_SOURCE"
else
    echo "ecbuild already cloned"
fi

if [ ! -d "$ECKIT_SOURCE" ]; then
    echo "Cloning eckit..."
    git clone --depth 1 --branch master https://github.com/ecmwf/eckit.git "$ECKIT_SOURCE"
else
    echo "eckit already cloned"
fi

if [ ! -d "$FCKIT_SOURCE" ]; then
    echo "Cloning fckit..."
    git clone --depth 1 --branch master https://github.com/ecmwf/fckit.git "$FCKIT_SOURCE"
else
    echo "fckit already cloned"
fi

if [ ! -d "$ATLAS_SOURCE/.git" ]; then
    if [ -e "$ATLAS_SOURCE" ]; then
        echo "Error: Atlas source exists but is not a Git checkout: $ATLAS_SOURCE"
        exit 1
    fi
    echo "Cloning Atlas $ATLAS_GIT_VERSION..."
    git clone --depth 1 --branch "$ATLAS_GIT_VERSION" \
        "$ATLAS_GIT_REPOSITORY" "$ATLAS_SOURCE"
else
    if [ -n "$(git -C "$ATLAS_SOURCE" status --porcelain)" ]; then
        echo "Error: Atlas source has local changes: $ATLAS_SOURCE"
        echo "Commit or remove those changes before selecting $ATLAS_GIT_VERSION."
        exit 1
    fi

    echo "Updating Atlas to $ATLAS_GIT_VERSION..."
    git -C "$ATLAS_SOURCE" fetch --depth 1 origin "$ATLAS_GIT_VERSION"
    if git -C "$ATLAS_SOURCE" show-ref --verify --quiet \
        "refs/heads/$ATLAS_GIT_VERSION"; then
        git -C "$ATLAS_SOURCE" checkout "$ATLAS_GIT_VERSION"
    else
        git -C "$ATLAS_SOURCE" checkout -b "$ATLAS_GIT_VERSION" FETCH_HEAD
    fi
    git -C "$ATLAS_SOURCE" merge --ff-only FETCH_HEAD
fi

# ========================================
# 2. Build ecbuild
# ========================================
echo ""
echo "========================================="
echo "Building ecbuild..."
echo "========================================="
mkdir -p "$ECBUILD_BUILD"
cd "$ECBUILD_BUILD"

if ! [ -f CMakeCache.txt ]; then
    cmake "$ECBUILD_SOURCE" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
        -DCMAKE_C_COMPILER="$CC" \
        -DCMAKE_CXX_COMPILER="$CXX" \
        2>&1 | tail -20
fi

echo "Building..."
make -j"$JOBS" 2>&1 | tail -5
echo "Installing..."
make install 2>&1 | tail -5

# ========================================
# 3. Build eckit
# ========================================
echo ""
echo "========================================="
echo "Building eckit..."
echo "========================================="
mkdir -p "$ECKIT_BUILD"
cd "$ECKIT_BUILD"

if ! [ -f CMakeCache.txt ]; then
    cmake "$ECKIT_SOURCE" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
        -DCMAKE_Fortran_COMPILER="$FC" \
        -DCMAKE_C_COMPILER="$CC" \
        -DCMAKE_CXX_COMPILER="$CXX" \
        -Decbuild_DIR="$INSTALL_PREFIX/lib/cmake/ecbuild" \
        -DCMAKE_PREFIX_PATH="$INSTALL_PREFIX" \
        -DENABLE_TESTS=OFF \
        -DENABLE_DOCS=OFF \
        -DENABLE_MPI=ON \
        2>&1 | tail -20
fi

echo "Building..."
make -j"$JOBS" 2>&1 | tail -5
echo "Installing..."
make install 2>&1 | tail -5

# ========================================
# 4. Build fckit (needs eckit)
# ========================================
echo ""
echo "========================================="
echo "Building fckit..."
echo "========================================="
mkdir -p "$FCKIT_BUILD"
cd "$FCKIT_BUILD"

if ! [ -f CMakeCache.txt ]; then
    cmake "$FCKIT_SOURCE" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
        -DCMAKE_Fortran_COMPILER="$FC" \
        -DCMAKE_C_COMPILER="$CC" \
        -DCMAKE_CXX_COMPILER="$CXX" \
    -Decbuild_DIR="$INSTALL_PREFIX/lib/cmake/ecbuild" \
    -DCMAKE_PREFIX_PATH="$INSTALL_PREFIX" \
        -Deckit_DIR="$INSTALL_PREFIX/lib/cmake/eckit" \
        -DENABLE_TESTS=OFF \
        -DENABLE_DOCS=OFF \
        2>&1 | tail -20
fi

echo "Building..."
make -j"$JOBS" 2>&1 | tail -5
echo "Installing..."
make install 2>&1 | tail -5

# ========================================
# 5. Build Atlas (needs eckit and fckit)
# ========================================
echo ""
echo "========================================="
echo "Building Atlas..."
echo "========================================="
mkdir -p "$ATLAS_BUILD"
cd "$ATLAS_BUILD"

if ! [ -f CMakeCache.txt ]; then
    cmake "$ATLAS_SOURCE" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
        -DCMAKE_Fortran_COMPILER="$FC" \
        -DCMAKE_C_COMPILER="$CC" \
        -DCMAKE_CXX_COMPILER="$CXX" \
    -Decbuild_DIR="$INSTALL_PREFIX/lib/cmake/ecbuild" \
    -DCMAKE_PREFIX_PATH="$INSTALL_PREFIX" \
        -Deckit_DIR="$INSTALL_PREFIX/lib/cmake/eckit" \
        -Dfckit_DIR="$INSTALL_PREFIX/lib/cmake/fckit" \
        -DENABLE_TESTS=OFF \
        -DENABLE_DOCS=OFF \
        -DENABLE_FORTRAN=ON \
        -DENABLE_OMP=OFF \
        -DENABLE_MPI=ON \
        2>&1 | tail -20
fi

echo "Building..."
make -j"$JOBS" 2>&1 | tail -10
echo "Installing..."
make install 2>&1 | tail -5

# ========================================
# 6. Build atlas-fesom plugin
# ========================================
echo ""
echo "========================================="
echo "Building atlas-fesom plugin..."
echo "========================================="

# Clone atlas-fesom if not already present
FESOM_PLUGIN_SOURCE="${DEPS_DIR}/atlas-fesom"
if [ ! -d "$FESOM_PLUGIN_SOURCE" ]; then
    echo "Cloning atlas-fesom from GitHub..."
    git clone https://github.com/ecmwf/atlas-fesom.git "$FESOM_PLUGIN_SOURCE" 2>&1 | tail -5
    cd atlas-fesom
    git checkout fix/ghost
else
    echo "atlas-fesom source already present"
fi

FESOM_PLUGIN_BUILD="${BUILD_BASE}/atlas-fesom"
mkdir -p "$FESOM_PLUGIN_BUILD"
cd "$FESOM_PLUGIN_BUILD"

if ! [ -f CMakeCache.txt ]; then
    cmake "$FESOM_PLUGIN_SOURCE" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$INSTALL_PREFIX" \
        -DCMAKE_Fortran_COMPILER="$FC" \
        -DCMAKE_C_COMPILER="$CC" \
        -DCMAKE_CXX_COMPILER="$CXX" \
        -DCMAKE_PREFIX_PATH="$INSTALL_PREFIX" \
        -Decbuild_DIR="$INSTALL_PREFIX/lib/cmake/ecbuild" \
        -Datlas_DIR="$INSTALL_PREFIX/lib/cmake/atlas" \
        -DENABLE_TESTS=OFF \
        -DENABLE_DOCS=OFF \
        2>&1 | tail -20
fi

echo "Building..."
make -j"$JOBS" 2>&1 | tail -10
echo "Installing..."
make install 2>&1 | tail -5

echo ""
echo "========================================="
echo "Build complete!"
echo "========================================="
echo "Install prefix: $INSTALL_PREFIX"
echo "Export this for use in configure.sh:"
echo "  export atlas_DIR=$INSTALL_PREFIX/lib/cmake/atlas"
echo ""
