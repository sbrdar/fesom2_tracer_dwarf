# Precision Support

The tracer dwarf supports three working precision modes controlled at compile time.

| | Double (dp) | Single (sp) | Half (hp) |
|---|:-:|:-:|:-:|
| `real` size | 8 bytes | 4 bytes | 2 bytes |
| Significant digits | ~15 | ~7 | ~3-4 |
| Range | ±1.8e308 | ±3.4e38 | ±6.5e4 |
| Epsilon | 2.2e-16 | 1.2e-7 | 9.8e-4 |
| GNU | yes | yes | no |
| Intel | yes | yes | no |
| NVIDIA | yes | yes | yes |
| Status | production | production | experimental |

## Building

```bash
# Double precision (default)
./configure.sh --compiler gnu --precision dp --clean --build

# Single precision
./configure.sh --compiler gnu --precision sp --clean --build

# Half precision (NVIDIA only)
./configure.sh --compiler nvidia --precision hp --clean --build
```

Or with CMake directly:

```bash
cmake .. -DUSE_SINGLE_PRECISION=ON    # single
cmake .. -DUSE_HALF_PRECISION=ON      # half (NVIDIA only)
cmake ..                               # double (default)
```

Build directories are named `build_<compiler>_<precision>`, allowing side-by-side builds.

## Implementation

### Working Precision (WP)

A three-way preprocessor in `lib/oce_modules.F90` sets the `WP` parameter:

```fortran
#ifdef USE_HALF_PRECISION
integer, parameter :: WP=2
#elif defined(USE_SINGLE_PRECISION)
integer, parameter :: WP=4
#else
integer, parameter :: WP=8
#endif
```

All real declarations throughout the codebase use `real(kind=WP)`. Literal constants use the `_WP` suffix (e.g., `1.e-3_WP`).

### Mesh Precision (MP)

A second parameter `MP = max(WP, 4)` ensures mesh and geometric arrays are always at least single precision. This prevents overflow in coordinate arithmetic (e.g., `b * r_earth` overflows FP16's ±65504 range).

- `T_MESH` arrays: `real(kind=MP)` -- coordinates, areas, volumes, geometry
- `T_TRACER_WORK` arrays: `real(kind=MP)` -- fluxes, FCT arrays, tendencies
- `T_TRACER_DATA` arrays: `real(kind=WP)` -- the actual tracer values under test
- Physical constants (`pi`, `r_earth`, `omega`): `real(kind=MP)`
- Diagnostic sums: wrapped with `dble()` to avoid accumulation overflow

When `WP >= 4` (SP or DP), `MP = WP` and all MP changes are no-ops.

### MPI Datatypes

`MPI_WP` maps to the correct MPI datatype per precision:

```fortran
#ifdef USE_SINGLE_PRECISION
integer, parameter :: MPI_WP = MPI_REAL
#else
integer, parameter :: MPI_WP = MPI_DOUBLE_PRECISION
#endif
```

For half precision, halo exchange calls are compiled out (`#if !defined(USE_HALF_PRECISION)`) since the analytic dwarf runs on a single rank.

### Compiler Flags

| Compiler | DP flag | SP flag | HP flag |
|----------|---------|---------|---------|
| GNU | `-fdefault-real-8` | (none, WP=4 via preprocessor) | not supported |
| Intel | `-r8` | `-r4` | not supported |
| NVIDIA | `-r8` | `-r4` | (none, WP=2 via preprocessor) |
| Cray | `-s real64` | `-s real32` | not supported |

### Source Files Modified

- `lib/oce_modules.F90` -- WP/MP parameters, MPI_WP
- `lib/gen_halo_exchange.F90` -- 58 hardcoded `real*8` replaced with `real(kind=WP)`; conditional exclusion of `real4`/`real8to4` gather procedures in SP mode
- `lib/oce_mesh.F90` -- `MPI_DOUBLE_PRECISION` replaced with `MPI_WP` (14 occurrences); local variables in mesh computation changed to MP
- `lib/gen_modules_partitioning.F90` -- `MPI_DOUBLE_PRECISION` replaced with `MPI_WP` (12 occurrences)
- `lib/MOD_MESH.F90` -- all arrays changed to `real(kind=MP)`
- `lib/MOD_TRACER.F90` -- work arrays to `real(kind=MP)`, data arrays stay `real(kind=WP)`
- `lib/gen_modules_rotate_grid.F90` -- added `trim_cyclic_mp` for MP-kind cyclic trimming
- `lib/hp_math_intrinsics.F90` -- SP-promoted wrappers for 10 intrinsics (HP only)

## Verification

### Test Results

**GNU DP/SP and NVIDIA DP/SP**: tracer sum perfectly conserved, 0 NaN at every step.

**NVIDIA HP**: all values finite, 0 NaN, sum drifts slightly (expected at FP16 precision).

### How to Verify Precision

**Check build output** -- CMake prints the precision mode:
```
-- Building with SINGLE precision (WP=4)
```

**Check runtime output** -- the executable prints WP at startup:
```
Working Precision (WP) =  4 bytes
*** SINGLE PRECISION MODE ***
```

**Compare binary size** -- SP shared library is ~20-30% smaller than DP.

**Compare memory usage**:
```bash
/usr/bin/time -v mpirun -np 1 ./bin/fesom_tracer_analytic 20 20 10 --periodic
# Look for "Maximum resident set size"
```

### Performance Notes

On small test meshes (e.g., pi mesh with 3140 nodes), runtime differences between SP and DP are small (5-15%) because the workload is I/O and MPI dominated. To see meaningful speedups:

- Use larger meshes (>100k nodes)
- Run more timesteps
- Measure specific advection kernels with `MPI_Wtime()`
- SP typically gives 1.1-1.5x CPU speedup, 1.5-2x GPU speedup at 50% memory

## Half Precision (FP16)

Half precision is experimental and only supported by NVIDIA's `nvfortran`.

### Limitations

- **Range**: ±65504 maximum. Ocean mesh areas (~1e14 m^2) and `r_earth` (6.4e6) overflow -- hence the MP design keeping mesh arrays at SP.
- **Precision**: ~3-4 digits. Temperature differences smaller than ~0.01 are lost.
- **Accumulation**: Summing many small values loses accuracy rapidly.
- **MPI**: No native MPI FP16 datatype. Halo exchange compiled out for npes=1.
- **Intrinsics**: `nvfortran` lacks FP16 overloads for `cos`, `sin`, `sqrt`, etc. The module `hp_math_intrinsics.F90` provides wrappers that promote to SP, compute, and demote back.
- **CPU performance**: No benefit (no hardware FP16 on x86 CPUs).
- **GPU performance**: 2-8x with tensor cores (Volta/Turing/Ampere+), minimal otherwise.

### Recommendation

Use double precision for production climate simulations. Use single precision for development and rapid testing. Half precision is only useful for GPU performance research and is not recommended for scientific results.
