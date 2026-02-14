# C and Python Interfaces

The tracer dwarf exposes a C-compatible interface (`iso_c_binding`) that can be called from Python via `ctypes`, or from any language with C FFI support. Built conditionally with `-DENABLE_TRACER_C_INTERFACE=ON`.

## Building

```bash
# Add to configure.sh or pass directly to CMake
./configure.sh --compiler gnu --precision dp --clean --build
# Then rebuild with C interface:
cd build_gnu_dp
cmake .. -DENABLE_TRACER_C_INTERFACE=ON
make -j$(nproc)

# Verify exported symbols
nm -D lib/libfesom_tracer_Fortran.so | grep "tracer_"
```

Source file: `lib/tracer_c_interface.F90`

## Function Reference

All functions return `int` status: 0 = success, non-zero = error. Error messages are written to stderr.

### MPI Lifecycle

| Function | C Signature | Description |
|----------|-------------|-------------|
| `tracer_init_mpi` | `int tracer_init_mpi(int mpi_comm)` | Initialize MPI (0 = MPI_COMM_WORLD) |
| `tracer_finalize_mpi` | `int tracer_finalize_mpi()` | Finalize MPI if we initialized it |

### Restart-Based Initialization

| Function | C Signature | Description |
|----------|-------------|-------------|
| `tracer_init` | `int tracer_init(char* restart_dir)` | Load mesh, partition, dynamics, tracers from binary restart files |
| `tracer_load_mesh_partition` | `int tracer_load_mesh_partition(char* restart_dir)` | Load only mesh and partition from restart files |
| `tracer_finalize` | `int tracer_finalize(char* output_dir)` | Write restart files and cleanup |

### Array-Based Initialization

These functions allow initializing the tracer system from NumPy arrays, without restart files. They must be called in the order listed.

| Function | C Signature | Description |
|----------|-------------|-------------|
| `tracer_set_mesh_dims` | `int tracer_set_mesh_dims(int nl, int nod2D, int edge2D, int elem2D)` | Set mesh dimensions |
| `tracer_set_levels` | `int tracer_set_levels(int* nlevels, int* ulevels, int nod2D)` | Set bottom/top level indices per node |
| `tracer_set_thickness` | `int tracer_set_thickness(double* hnode, int nl, int nod2D)` | Set layer thickness (nl-1 x nod2D) |
| `tracer_set_partition` | `int tracer_set_partition(int myDim_nod2D, int eDim_nod2D)` | Set partition dimensions |
| `tracer_allocate_tracers` | `int tracer_allocate_tracers(int num_tracers, int nl, int nod2D, int AB_order)` | Allocate tracer data and work arrays |
| `tracer_set_values` | `int tracer_set_values(int tracer_id, double* values, int nl, int nod2D)` | Set tracer field (nl-1 x nod2D, 1-based ID) |
| `tracer_set_velocity` | `int tracer_set_velocity(double* uv, double* w, int nl, int nod2D, int edge2D)` | Set velocity from arrays |
| `tracer_set_velocity_zero` | `int tracer_set_velocity_zero(int nl, int nod2D, int edge2D)` | Allocate velocity arrays, set to zero |
| `tracer_init_complete` | `int tracer_init_complete()` | Mark initialization done |

### Runtime

| Function | C Signature | Description |
|----------|-------------|-------------|
| `tracer_advect_step` | `int tracer_advect_step(double dt)` | Perform one advection time step |
| `tracer_get_stats` | `int tracer_get_stats(double* tmin, double* tmax, double* tsum)` | Get tracer min, max, sum |
| `tracer_get_size` | `int tracer_get_size(int* nz, int* nn)` | Get tracer array dimensions |
| `tracer_get_values` | `int tracer_get_values(double* values, int nz, int nn)` | Copy tracer array to caller |
| `tracer_cleanup` | `int tracer_cleanup()` | Reset module state |
| `tracer_run_workflow` | `int tracer_run_workflow(char* restart_dir, int nsteps, double dt)` | Complete init-run-finalize in one call |

## Initialization Order

### Restart-based (simplest)

```
tracer_init_mpi(0)
tracer_init("path/to/restart")       # loads everything
tracer_advect_step(dt)               # repeat as needed
tracer_finalize("path/to/output")
tracer_finalize_mpi()
```

### Array-based (no restart files)

```
tracer_init_mpi(0)
tracer_set_mesh_dims(nl, nod2D, edge2D, elem2D)
tracer_set_levels(nlevels, ulevels, nod2D)
tracer_set_thickness(hnode, nl, nod2D)
tracer_set_partition(myDim_nod2D, eDim_nod2D)
tracer_allocate_tracers(num_tracers, nl, nod2D, AB_order)
tracer_set_values(tracer_id, values, nl, nod2D)   # per tracer
tracer_set_velocity_zero(nl, nod2D, edge2D)        # or tracer_set_velocity()
tracer_init_complete()
tracer_advect_step(dt)                              # repeat as needed
tracer_finalize_mpi()
```

### Hybrid (load mesh from files, set tracers from arrays)

```
tracer_init_mpi(0)
tracer_load_mesh_partition("path/to/restart")
tracer_allocate_tracers(num_tracers, nl, nod2D, AB_order)
tracer_set_values(tracer_id, values, nl, nod2D)
tracer_init_complete()
tracer_advect_step(dt)
tracer_finalize_mpi()
```

## Array Layouts

All arrays use **Fortran column-major order**. In NumPy, use `order='F'` and `dtype=np.float64`.

| Array | Shape | Notes |
|-------|-------|-------|
| Tracer values | `(nl-1, nod2D)` | Level x node |
| Horizontal velocity `uv` | `(nl-1, edge2D, 2)` | Level x edge x component (u,v) |
| Vertical velocity `w` | `(nl, nod2D)` | Level x node |
| Layer thickness `hnode` | `(nl-1, nod2D)` | Level x node |
| Level indices `nlevels`, `ulevels` | `(nod2D,)` | Integer, per node |

## Python Usage

### Loading the Library

```python
import ctypes
import numpy as np

lib = ctypes.CDLL("build_gnu_dp/lib/libfesom_tracer_Fortran.so")
```

### Python Wrapper Class

The `python/demo_tracer.py` file provides a `FESOMTracerAdvection` class:

```python
from demo_tracer import FESOMTracerAdvection

tracer = FESOMTracerAdvection()
tracer.initialize_mpi(0)
tracer.initialize("../fesom_bin_restart")

for step in range(10):
    tracer.advect_step(1.0e-3)
    stats = tracer.get_statistics()
    print(f"Step {step}: min={stats['min']:.4f}, max={stats['max']:.4f}, sum={stats['sum']:.2e}")

tracer.finalize("../fesom_bin_restart")
tracer.finalize_mpi()
```

### Array Initialization Example

```python
tracer = FESOMTracerAdvection()
tracer.initialize_mpi(0)

nl, nod2D, edge2D, elem2D = 10, 100, 200, 90

# Mesh setup
tracer.lib.tracer_set_mesh_dims(nl, nod2D, edge2D, elem2D)

nlevels = np.full(nod2D, nl - 1, dtype=np.int32)
ulevels = np.ones(nod2D, dtype=np.int32)
tracer.lib.tracer_set_levels(nlevels, ulevels, nod2D)

thickness = np.ones((nl - 1, nod2D), dtype=np.float64, order='F') * 100.0
tracer.lib.tracer_set_thickness(thickness, nl, nod2D)

tracer.lib.tracer_set_partition(myDim_nod2D=nod2D, eDim_nod2D=0)

# Tracers
tracer.lib.tracer_allocate_tracers(num_tracers=2, nl=nl, nod2D=nod2D, AB_order=2)

temperature = np.zeros((nl - 1, nod2D), dtype=np.float64, order='F')
for k in range(nl - 1):
    temperature[k, :] = 20.0 - 15.0 * (k / (nl - 2))
tracer.lib.tracer_set_values(tracer_id=1, values=temperature, nl=nl, nod2D=nod2D)

salinity = np.full((nl - 1, nod2D), 35.0, dtype=np.float64, order='F')
tracer.lib.tracer_set_values(tracer_id=2, values=salinity, nl=nl, nod2D=nod2D)

# Velocity
tracer.lib.tracer_set_velocity_zero(nl, nod2D, edge2D)
tracer.lib.tracer_init_complete()

tracer.finalize_mpi()
```

## Design Notes

- **Module-level storage**: The Fortran module holds one global instance of each derived type (`t_mesh`, `t_partit`, `t_tracer`, `t_dyn`). Only one tracer system can be active at a time.
- **Memory ownership**: Fortran allocates and owns all arrays. Python passes data that gets copied into Fortran arrays.
- **MPI**: The interface handles MPI init/finalize. Works with `mpirun -np N` for parallel execution.
- **Type conversion**: C `double` maps to Fortran `real(c_double)`. Automatic conversion to/from `real(kind=WP)` internally.

## Known Limitations

- **Advection requires mesh connectivity**: Pure array initialization sets up basic arrays but advection needs edge/element connectivity. Use `tracer_load_mesh_partition()` or full `tracer_init()` for advection runs.
- **Single instance**: Module-level storage limits to one active tracer system.
- **Statistics**: `tracer_get_stats` returns stats for tracer 1 only. Per-tracer stats would need `tracer_get_stats_by_id()`.
- **No native MPI FP16**: Half precision builds compile out halo exchange (single-rank only).

## Files

| File | Description |
|------|-------------|
| `lib/tracer_c_interface.F90` | C interface implementation |
| `python/demo_tracer.py` | Python wrapper class (`FESOMTracerAdvection`) |
| `python/test_python_init.py` | Array initialization test script |
| `python/verify_c_interface.py` | Symbol verification script |
