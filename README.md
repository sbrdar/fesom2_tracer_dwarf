# FESOM2 Tracer Advection Dwarf

A standalone tracer advection kernel extracted from FESOM2. It tests the horizontal and vertical tracer transport (MUSCL, FCT, upwind, PPM schemes) on unstructured triangular meshes with multi-precision support (double, single, half precision).

## Quick Start

Build and run the analytic dwarf (no input files needed):

```bash
./configure.sh --compiler gnu --precision dp --clean --build
cd build_gnu_dp
./run.sh 1 20 20 10 --periodic
```

This generates a 20x20 doubly-periodic triangular mesh with 10 vertical levels and runs 10 advection steps on a single MPI rank.

Build the Atlas dependencies, then configure the dwarf with Atlas enabled:

```bash
./build_atlas.sh --compiler gnu
./configure.sh --compiler gnu --precision dp --clean --build --atlas
```

`build_atlas.sh` delegates to Atlas's `tools/install.sh` with all dependencies,
Fortran, LZ4, and atlas-fesom enabled. The atlas-fesom dependency path also
installs METIS. If the Atlas source is absent, the script clones it into
`atlas_deps/atlas`. With `--atlas`, `configure.sh` enables atlas-fesom caching in
the generated run wrapper, allowing automatic download of requested FESOM grids.
The resulting binaries use the standard non-Atlas path by default.
Set `ATLAS_FESOM=1` at runtime to use Atlas:

```bash
cd build_gnu_dp_atlas
ATLAS_FESOM=1 ATLAS_GRID=fesom-pi ATLAS_USE_FESOM_DIST=0 mpirun -np 8 ./bin/fesom_tracer_mesh_init
```

`ATLAS_GRID` selects the Atlas grid and defaults to `fesom-pi` when unset or
empty.
`ATLAS_USE_FESOM_DIST=1` reads in the original FESOM paritioning from a file and applies it to Atlas generated grid, if `ATLAS_FESOM=1` is used.

The mesh-file driver also accepts Atlas structured grids on one MPI rank:

```bash
ATLAS_FESOM=1 ATLAS_GRID=O64 ./bin/fesom_tracer_mesh_init
```

## Running

### Analytic Dwarf (recommended for testing)

Generates a rectangular triangular mesh in-memory. No mesh files or restart data required.

```bash
# Usage: ./run.sh NP nx ny nl [options]
#   NP  = number of MPI processes (always 1 for analytic)
#   nx  = grid points in x (default: 20)
#   ny  = grid points in y (default: 20)
#   nl  = vertical levels (default: 10, minimum: 3)

# Basic run
./run.sh 1 20 20 10 --periodic

# Larger mesh with output
./run.sh 1 50 50 20 --periodic --save-scalars

# Save mesh arrays for visualization
./run.sh 1 40 40 10 --periodic --save-mesh --save-scalars

# Direct invocation (equivalent)
mpirun -np 1 ./bin/fesom_tracer_analytic 40 40 10 --periodic --save-scalars
```

Options:
- `--periodic` : Doubly-periodic boundary conditions (default: closed-wall)
- `--save-mesh` : Write mesh arrays as raw binary for Python visualization
- `--save-scalars` : Write tracer fields at each step
- `--output-dir DIR` : Output directory (default: `mesh_output`)

Expected output shows perfectly conserved tracer sum at each step:
```
Step    1: T min, max, sum =   0.246394E-01  0.199736E+02  0.3249000000E+05  NaN:     0
Step    2: T min, max, sum =   0.322013E-01  0.199642E+02  0.3249000000E+05  NaN:     0
...
```

### Mesh-Init Dwarf

Reads global mesh files (not restart files) and sets custom tracer values before
advection. A single-rank run derives identity ownership from the global mesh;
multi-rank runs require matching `dist_N/` partition files.

```bash
mpirun -np 1 ./bin/fesom_tracer_mesh_init
```

### Restart Dwarf

Reads binary restart files (mesh, partition, dynamics, tracers), runs advection, and writes updated restarts. Requires restart data in `../fesom_bin_restart/np<N>/`.

```bash
mpirun -np <N> ./bin/fesom_tracer
```

## Building

### Using `configure.sh`

```bash
./configure.sh [OPTIONS]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `--compiler` | `gnu`, `intel`, `nvidia` | `gnu` | Fortran compiler |
| `--precision` | `dp`, `sp`, `hp` | `dp` | Working precision |
| `--openacc` | | off | Enable OpenACC (NVIDIA only) |
| `--build-type` | `Release`, `Debug` | `Release` | CMake build type |
| `--clean` | | | Remove build dir before configuring |
| `--build` | | | Run `make` after configuring |

Build directories are named `build_<compiler>_<precision>` (e.g., `build_gnu_dp`, `build_nvidia_sp`).

Examples:

```bash
# GNU double precision (default), configure + build
./configure.sh --clean --build

# GNU single precision
./configure.sh --compiler gnu --precision sp --clean --build

# Intel double precision
./configure.sh --compiler intel --precision dp --clean --build

# NVIDIA with OpenACC
./configure.sh --compiler nvidia --openacc --clean --build

# NVIDIA half precision (experimental, requires nvfortran)
./configure.sh --compiler nvidia --precision hp --clean --build
```

### Visualization

Run the analytic dwarf with `--save-mesh` and `--save-scalars` to write binary output, then use the plotting scripts in `visualization/` to inspect results.

```bash
# Run and save output
cd build_gnu_dp
./run.sh 1 30 30 10 --periodic --save-mesh --save-scalars

# Plot the mesh (triangulation, bathymetry, vertical levels)
python ../visualization/plot_mesh.py mesh_output

# Plot tracer fields at each timestep (surface temperature/salinity)
python ../visualization/plot_scalars.py mesh_output

# Compare two precision runs side by side
cd ..
./configure.sh --compiler gnu --precision sp --clean --build
cd build_gnu_sp
./run.sh 1 30 30 10 --periodic --save-mesh --save-scalars --output-dir mesh_output
cd ..
python visualization/plot_compare_precision.py --diff \
    build_gnu_dp/mesh_output build_gnu_sp/mesh_output
```

The scripts require `numpy` and `matplotlib`. Output is saved as PNG images (and animated GIFs for `plot_scalars.py`) in the mesh output directory.

### Manual CMake

```bash
mkdir build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_C_COMPILER=gcc
make -j$(nproc)
```

CMake options:
- `-DUSE_SINGLE_PRECISION=ON` : Single precision (WP=4)
- `-DUSE_HALF_PRECISION=ON` : Half precision (WP=2, NVIDIA only)
- `-DENABLE_OPENACC=ON` : OpenACC GPU offloading
- `-DENABLE_TRACER_C_INTERFACE=ON` : Build C interface for Python binding

### Compiler/Precision Support

| Compiler | Double (dp) | Single (sp) | Half (hp) |
|----------|:-----------:|:-----------:|:---------:|
| GNU (gfortran) | yes | yes | no |
| Intel (ifx) | yes | yes | no |
| NVIDIA (nvfortran) | yes | yes | yes |

Half precision is experimental. NVIDIA's `nvfortran` supports `real(kind=2)` but lacks FP16 intrinsic overloads; the dwarf includes `hp_math_intrinsics.F90` to provide SP-promoted wrappers.

### Build Products

```
build_<compiler>_<precision>/
  bin/
    fesom_tracer            # restart-based driver
    fesom_tracer_analytic   # analytic mesh driver
    fesom_tracer_mesh_init  # mesh-file driver
  lib/
    libfesom_tracer_Fortran.so  # shared library
```

## Further Documentation

- [Architecture](docs/ARCHITECTURE.md) -- directory layout and source file descriptions
- [Precision](docs/PRECISION.md) -- mixed precision design (WP/MP), compiler support, FP16 details
- [C and Python Interfaces](docs/INTERFACES.md) -- C-compatible API, Python wrapper, array layouts
- [Visualization](docs/VISUALIZATION.md) -- plotting scripts setup, usage, and multi-precision comparison

## Note for AI editors

See [docs/CLAUDE.md](docs/CLAUDE.md) for codebase conventions, key file locations, and precision pitfalls relevant to AI-assisted editing.
