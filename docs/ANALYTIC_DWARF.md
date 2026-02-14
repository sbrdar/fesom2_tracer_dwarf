# FESOM2 Analytic Tracer Dwarf

A self-contained tracer advection test that generates a rectangular triangular mesh entirely in-memory, with no external mesh files required. Supports both closed-wall and doubly-periodic boundary conditions.

## Building

### Quick Start (GNU, double precision)

```bash
./configure.sh --build
```

### Using `configure.sh`

The `configure.sh` script handles compiler detection, MPI setup, and precision flags. Build directories are named `build_<compiler>_<precision>`.

```bash
./configure.sh [OPTIONS]
```

| Option | Values | Default | Description |
|--------|--------|---------|-------------|
| `--compiler` | `gnu`, `intel`, `nvidia` | `gnu` | Fortran compiler family |
| `--precision` | `dp`, `sp`, `hp` | `dp` | Working precision (double, single, half) |
| `--build` | - | off | Also run `make` after configuring |
| `--clean` | - | off | Remove build directory before configuring |
| `--openacc` | - | off | Enable OpenACC (NVIDIA only) |
| `--build-type` | `Release`, `Debug` | `Release` | CMake build type |

### Build Examples

```bash
# GNU double precision (default)
./configure.sh --build                                    # -> build_gnu_dp/

# GNU single precision
./configure.sh --compiler gnu --precision sp --build      # -> build_gnu_sp/

# Intel double precision
./configure.sh --compiler intel --precision dp --build    # -> build_intel_dp/

# Intel single precision
./configure.sh --compiler intel --precision sp --build    # -> build_intel_sp/

# NVIDIA double precision
./configure.sh --compiler nvidia --precision dp --build   # -> build_nvidia_dp/

# NVIDIA half precision
./configure.sh --compiler nvidia --precision hp --build   # -> build_nvidia_hp/
```

A `build` symlink is created pointing to the latest build directory.

### Compiler Details

| Compiler | ID in CMake | Precision Support | Real-size Flag | Notes |
|----------|-------------|-------------------|----------------|-------|
| GNU (gfortran) | `GNU` | DP, SP | `-fdefault-real-8` (DP only) | No `real(kind=2)` support |
| Intel (ifx) | `IntelLLVM` | DP, SP | `-r8` / `-r4` | No `real(kind=2)` support |
| NVIDIA (nvfortran) | `NVHPC` | DP, SP, HP | `-r8` / `-r4` / (none for HP) | HP uses `real(kind=2)` via preprocessor |

### Manual CMake Build

```bash
mkdir -p build && cd build
cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_Fortran_COMPILER=gfortran \
    -DCMAKE_C_COMPILER=gcc \
    -DUSE_SINGLE_PRECISION=OFF \
    -DUSE_HALF_PRECISION=OFF \
    -DENABLE_OPENACC=OFF
make -j$(nproc)
```

### Precision CMake Options

| Option | Effect |
|--------|--------|
| Neither `USE_SINGLE_PRECISION` nor `USE_HALF_PRECISION` | Double precision (WP=8, default) |
| `-DUSE_SINGLE_PRECISION=ON` | Single precision (WP=4) |
| `-DUSE_HALF_PRECISION=ON` | Half precision (WP=2, NVIDIA only) |

Both options cannot be enabled simultaneously.

## Usage

Each build produces its executable at `build_<compiler>_<precision>/bin/fesom_tracer_analytic`. Run from within the build directory so the shared library resolves correctly via rpath.

```bash
cd build_gnu_dp
mpirun -np 1 ./bin/fesom_tracer_analytic [nx] [ny] [nl] [OPTIONS]
```

### Positional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `nx`     | 20      | Grid points in x direction (must be >= 2) |
| `ny`     | 20      | Grid points in y direction (must be >= 2) |
| `nl`     | 10      | Number of vertical levels (must be >= 3) |

### Options

| Flag             | Description |
|------------------|-------------|
| `--save-mesh`    | Write mesh geometry as binary files |
| `--save-scalars` | Write tracer fields at every timestep |
| `--periodic`     | Use doubly-periodic boundary conditions (default: closed-wall) |
| `--output-dir DIR` | Output directory for mesh/scalar files (default: `mesh_output`) |

Flags can appear in any position and in any order.

### Fixed Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Domain Lx | 100 km | Domain size in x |
| Domain Ly | 100 km | Domain size in y |
| Max depth | 1000 m | Uniform flat bottom |
| Timesteps | 10 | Number of advection steps |
| dt | 1000 s | Timestep size |
| Velocity | (0.1, 0.0) m/s | Uniform eastward current |
| Advection | UPW1 | First-order upwind, horizontal + vertical |

### Examples

```bash
# Default 20x20 grid, closed-wall boundaries
cd build_gnu_dp
mpirun -np 1 ./bin/fesom_tracer_analytic

# Larger mesh with periodic boundaries and output
mpirun -np 1 ./bin/fesom_tracer_analytic 50 50 10 --save-mesh --save-scalars --periodic

# Save output to a named directory (useful for comparing builds)
mpirun -np 1 ./bin/fesom_tracer_analytic 50 50 10 --save-scalars --periodic --output-dir /tmp/output_gnu_dp

cd ../build_gnu_sp
mpirun -np 1 ./bin/fesom_tracer_analytic 50 50 10 --save-scalars --periodic --output-dir /tmp/output_gnu_sp
```

## Multi-Precision Comparison

Use `--output-dir` to save runs from different builds to separate directories, then compare with `plot_scalars.py --compare`:

```bash
# Run DP and SP builds with same parameters
cd build_gnu_dp && mpirun -np 1 ./bin/fesom_tracer_analytic 50 50 10 --save-mesh --save-scalars --periodic --output-dir /tmp/cmp_dp
cd ../build_gnu_sp && mpirun -np 1 ./bin/fesom_tracer_analytic 50 50 10 --save-mesh --save-scalars --periodic --output-dir /tmp/cmp_sp

# Compare side-by-side with difference panel
python plot_scalars.py --compare /tmp/cmp_dp /tmp/cmp_sp
```

The comparison mode produces 3-panel plots (field 1, field 2, difference) for each tracer at each timestep, along with per-step difference statistics (max|diff|, mean|diff|, RMS).

## Mesh Topology

The domain is a Cartesian rectangle `[0, Lx] x [0, Ly]` triangulated into right triangles. Each grid cell `(i,j)` is split into two triangles:

```
(i,j+1)-----(i+1,j+1)
  | \  upper  |
  |   \       |
  |     \     |
  |  lower \  |
  |         \ |
(i,j)------(i+1,j)
```

### Non-Periodic (Closed-Wall) Mode

- **Nodes:** `nx * ny`
- **Elements:** `2 * (nx-1) * (ny-1)`
- **Edges:** `(nx-1)*ny + nx*(ny-1) + (nx-1)*(ny-1)` total, of which `2*(nx-1) + 2*(ny-1)` are boundary edges

Boundary edges have `edge_tri(2,n) == 0`, which zeros out the flux contribution from the missing side. This implements a no-normal-flux (closed-wall) boundary condition.

### Periodic Mode

- **Real nodes:** `(nx-1) * (ny-1)` (grid points i=1..nx-1, j=1..ny-1)
- **Elements:** `2 * (nx-1) * (ny-1)`
- **Edges:** `3 * (nx-1) * (ny-1)`, all internal (zero boundary edges)

The i=nx column wraps to i=1, and the j=ny row wraps to j=1.

## Tracer Initialization

Two tracers are initialized:

**Temperature** (tracer 1): Sinusoidal 4-lobe pattern with vertical stratification
```
T(x, y, z) = 15 + 5 * sin(2*pi*x/Lx) * cos(2*pi*y/Ly) - 10*z_norm
```

**Salinity** (tracer 2): Gaussian blob at domain center with depth decay
```
S(x, y, z) = 34 + 2 * exp(-r^2 / (0.1*Lx)^2) * (1 - z_norm)
```

where `r = sqrt((x - Lx/2)^2 + (y - Ly/2)^2)` and `z_norm` goes from 0 (surface) to 1 (bottom).

## Periodic Boundary Implementation

The doubly-periodic boundary condition is implemented using a temporary halo (ghost) node approach.

### Halo Node Strategy

During mesh construction, temporary halo nodes are created with shifted coordinates so that FESOM's existing mesh geometry routines (`mesh_areas`, `mesh_auxiliary_arrays`) compute correct areas and edge properties for wrap-around elements:

| Node type | Count | Coordinate shift |
|-----------|-------|-------------------|
| Real nodes | `(nx-1)*(ny-1)` | None |
| Right halos | `ny-1` | Copies of i=1 column with x += Lx |
| Top halos | `nx-1` | Copies of j=1 row with y += Ly |
| Corner halo | 1 | Copy of (1,1) with x += Lx, y += Ly |

All edges are internal -- every edge has two adjacent elements due to periodic wrapping. This means zero boundary edges, and the advection code accumulates flux from both sides of every edge.

### Post-Remap

After `mesh_areas` and `mesh_auxiliary_arrays` run on the expanded mesh (real + halo nodes), a post-remap step:

1. **Fixes `edge_cross_dxdy`** using minimum-image convention: for wrap-around edges, the vector from edge center to element centroid spans the entire domain instead of wrapping correctly. If any component exceeds half the domain size, the full domain size is subtracted.

2. **Merges halo node areas** into real nodes: `area(real) += area(halo)`, `areasvol(real) += areasvol(halo)`.

3. **Recomputes derived quantities**: `area_inv`, `areasvol_inv`, `mesh_resolution`.

4. **Remaps connectivity arrays**: replaces halo node indices with real node indices in `elem2D_nodes` and `edges`.

5. **Sets `eDim_nod2D = 0`**: subsequent allocations (tracers, hnode, etc.) use only real nodes.

After remap, no periodic halo exchange is needed at runtime -- all flux accumulation happens at the correct real node indices with correct total areas.

## Output Files

### Mesh Output (`--save-mesh`)

Written to the output directory (default `mesh_output/`, or `--output-dir DIR`):

| File | Format | Shape | Description |
|------|--------|-------|-------------|
| `mesh_info.txt` | Text | - | `nod2D`, `elem2D`, `nl` |
| `coord_nod2D.bin` | float64, Fortran order | (2, nod2D) | Node x,y coordinates in meters |
| `elem2D_nodes.bin` | int32, Fortran order | (3, elem2D) | Triangle connectivity, 0-based |
| `zbar.bin` | float64 | (nl,) | Level interface depths |
| `depth.bin` | float64 | (nod2D,) | Bottom depth at each node |

### Scalar Output (`--save-scalars`)

Written to the output directory:

| File | Format | Shape | Description |
|------|--------|-------|-------------|
| `scalar_info.txt` | Text | - | `nod2D`, `nz`, `nsteps`, `num_tracers`, `wp`, tracer names |
| `temperature.bin` | float64, Fortran order | (nsteps+1) frames of (nz, nod2D) | Temperature at each step |
| `salinity.bin` | float64, Fortran order | (nsteps+1) frames of (nz, nod2D) | Salinity at each step |

Each binary file contains `nsteps + 1` frames (initial condition + one per timestep). Each frame is a Fortran column-major `(nz, nod2D)` array written sequentially. Binary output is always float64 regardless of working precision.

The `wp` field in `scalar_info.txt` records the working precision in bytes (8=DP, 4=SP, 2=HP), used by `plot_scalars.py` for labeling.

## Plotting Utilities

### `plot_mesh.py`

Plots the surface triangulation.

```bash
python plot_mesh.py [mesh_output_dir]
```

Default directory: `mesh_output`. Saves `mesh_plot.png` in the output directory.

### `plot_scalars.py`

Plots surface layer (k=0) of each tracer at each timestep and creates animated GIFs.

**Single directory mode:**

```bash
python plot_scalars.py [mesh_output_dir]
```

**Comparison mode** (side-by-side with difference panel):

```bash
python plot_scalars.py --compare dir1 dir2 [--output-dir compare_output/]
```

Comparison mode produces:
- 3-panel figures per step: field from dir1, field from dir2, difference (dir1 - dir2)
- Per-step difference statistics printed to console: max|diff|, mean|diff|, RMS
- Animated GIFs of the comparison
- Precision labels (DP/SP/HP) read from `scalar_info.txt`

**Output:**
- Per-step PNG images in `<dir>/scalar_plots/` (single mode) or `<output_dir>/compare_plots/` (compare mode)
- Animated GIFs (requires Pillow)

**Periodic mesh handling:** Wrap-around triangles (those connecting nodes at opposite edges of the domain) are automatically detected and masked. A triangle is masked if any edge span exceeds half the domain extent. This prevents rendering artifacts from large diagonal triangles.

### Dependencies

- `numpy` (required)
- `matplotlib` (required)
- `Pillow` (optional, for GIF animation)

## Half Precision (FP16) Support

Half precision is supported only with NVIDIA's nvfortran compiler (`--compiler nvidia --precision hp`).

### Compiler Support

| Compiler | Half Precision | Reason |
|----------|---------------|--------|
| GNU (gfortran) | Not supported | `real(kind=2)` is not a valid kind |
| Intel (ifx) | Not supported | `real(kind=2)` is not a valid kind; `-r2` flag unknown |
| NVIDIA (nvfortran) | Builds successfully | Supports `real(kind=2)` natively |

### Intrinsic Math Function Wrappers

nvfortran supports `real(kind=2)` for storage and arithmetic, but most intrinsic math functions (`cos`, `sin`, `tan`, `sqrt`, `abs`, `asin`, `atan2`, `sign`, `exp`, `log`) lack FP16 overloads. The module `src/hp_math_intrinsics.F90` provides elemental wrappers that:

1. Promote `real(2)` arguments to `real(4)` (single precision)
2. Call the intrinsic
3. Demote the result back to `real(2)`

The wrappers use Fortran generic interfaces with the same names as the intrinsics, so existing code only needs `use hp_math_intrinsics` (guarded by `#ifdef USE_HALF_PRECISION`) to get FP16 support transparently.

Source files with the conditional `use`:
- `src/gen_modules_rotate_grid.F90`
- `src/oce_mesh.F90`
- `src/oce_muscl_adv.F90`
- `src/analytic_mesh.F90`
- `src/oce_adv_tra_hor.F90` (3 subroutines)
- `src/oce_adv_tra_ver.F90` (5 subroutines)
- `dwarf_ini/fesom_analytic.F90`

### Runtime Limitation

FP16 (IEEE 754 half precision) has a maximum representable value of ~65,504 and minimum positive normal of ~6.1e-5. Physical constants in this code exceed FP16 range:

| Constant | Value | FP16 Status |
|----------|-------|-------------|
| `r_earth` | 6,367,500 m | Overflows to Inf |
| `Lx` | 100,000 m | Overflows to Inf |
| `surf_relax_S` denominator | 5,184,000 | Overflows to Inf |

This causes all downstream values (coordinates, areas, tracers) to become NaN. Making HP produce correct results at runtime would require fully non-dimensionalized physics, which is beyond the current scope. The HP build demonstrates compiler toolchain support and the intrinsic wrapper infrastructure.

### Additional Code Changes for HP

- `src/oce_adv_tra_fct.F90`: Replaced legacy Fortran 77 `dmax1`/`dmin1` intrinsics with standard `max`/`min` (these are DP-specific and have no FP16 overload).

## Shared Library and rpath

Each build produces a shared library `lib/libfesom_tracer_Fortran.so`. The executable uses `$ORIGIN/../lib` rpath to locate it. This means:

- Always run the executable from its own build directory
- Do NOT copy executables between build directories (precision mismatch between executable and library causes segfaults)
- Use `--output-dir` to redirect output to separate directories for comparison

## Files Modified

| File | Changes |
|------|---------|
| `CMakeLists.txt` | `apply_fesom_compile_flags()` function for unified compiler/precision flag handling; GNU/Intel/IntelLLVM/NVHPC/Cray/PGI support; HP precision option |
| `configure.sh` | Multi-compiler build script with `--compiler`, `--precision`, `--build`, `--clean`, `--openacc` |
| `src/hp_math_intrinsics.F90` | New: FP16-to-SP elemental wrappers for 10 intrinsic math functions (guarded by `#ifdef USE_HALF_PRECISION`) |
| `src/analytic_mesh.F90` | Added `periodic` parameter; periodic node/element/edge generation with halos; `pnode()` helper; `periodic_post_remap()` subroutine; conditional `use hp_math_intrinsics` |
| `src/oce_adv_tra_fct.F90` | `dmax1`/`dmin1` replaced with standard `max`/`min` |
| `src/gen_modules_rotate_grid.F90` | Conditional `use hp_math_intrinsics` |
| `src/oce_mesh.F90` | Conditional `use hp_math_intrinsics` |
| `src/oce_muscl_adv.F90` | Conditional `use hp_math_intrinsics` |
| `src/oce_adv_tra_hor.F90` | Conditional `use hp_math_intrinsics` (3 subroutines) |
| `src/oce_adv_tra_ver.F90` | Conditional `use hp_math_intrinsics` (5 subroutines) |
| `src/mesh_output.F90` | Added `wp=` to `scalar_info.txt` |
| `dwarf_ini/fesom_analytic.F90` | Added `--periodic`, `--output-dir` CLI flags; conditional `use hp_math_intrinsics` |
| `plot_scalars.py` | Added `--compare` mode, `--output-dir`, precision labels (DP/SP/HP), wrap-around triangle masking |

## Verification

```bash
# GNU DP: periodic, constant T sum
cd build_gnu_dp
mpirun -np 1 ./bin/fesom_tracer_analytic 50 50 10 --save-scalars --periodic --output-dir /tmp/v_dp

# GNU SP: periodic, T sum shows SP rounding
cd ../build_gnu_sp
mpirun -np 1 ./bin/fesom_tracer_analytic 50 50 10 --save-scalars --periodic --output-dir /tmp/v_sp

# Compare DP vs SP
python plot_scalars.py --compare /tmp/v_dp /tmp/v_sp

# Cross-compiler comparison (DP builds differ at ~1e-15)
cd ../build_intel_dp
mpirun -np 1 ./bin/fesom_tracer_analytic 50 50 10 --save-scalars --periodic --output-dir /tmp/v_intel_dp
python plot_scalars.py --compare /tmp/v_dp /tmp/v_intel_dp
```

Expected behavior:
- **DP builds**: T sum constant to machine epsilon (~1e-15); cross-compiler differences at ~1e-15 level
- **SP builds**: T sum shows small fluctuations from SP rounding; differences vs DP at ~1e-7 level
- **HP builds**: All NaN due to FP16 overflow of physical constants (build succeeds, runtime limited by FP16 range)
