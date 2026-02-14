# Architecture

```
fesom2_dwarf_tracer/
  CMakeLists.txt          # top-level build configuration
  configure.sh            # compiler/precision selection wrapper
  README.md

  src/                    # executable entry points
    fesom_analytic.F90        # analytic mesh driver (no I/O)
    fesom_mesh_init.F90       # mesh-file driver
    fesom.F90                 # restart-file driver

  lib/                    # shared library sources (libfesom_tracer_Fortran.so)
    oce_modules.F90           # WP/MP precision parameters, physical constants
    hp_math_intrinsics.F90    # FP16 math wrappers (half precision only)
    gen_modules_config.F90    # model configuration / namelists

    MOD_MESH.F90              # mesh derived type (T_MESH)
    MOD_PARTIT.F90            # partition derived type (T_PARTIT)
    MOD_TRACER.F90            # tracer derived types (T_TRACER, T_TRACER_WORK, T_TRACER_DATA)
    MOD_DYN.F90               # dynamics derived type (T_DYN)

    oce_adv_tra_driver.F90    # top-level advection driver (do_oce_adv_tra)
    oce_adv_tra_hor.F90       # horizontal advection (MUSCL, MFCT, UPW1)
    oce_adv_tra_ver.F90       # vertical advection (QR4C, CDIFF, PPM, UPW1)
    oce_adv_tra_fct.F90       # FCT limiter
    oce_muscl_adv.F90         # MUSCL gradient reconstruction

    oce_mesh.F90              # mesh geometry computation
    analytic_mesh.F90         # analytic mesh generation
    mesh_output.F90           # binary mesh/scalar output

    gen_modules_rotate_grid.F90   # coordinate rotation utilities
    gen_modules_partitioning.F90  # mesh partitioning
    gen_halo_exchange.F90         # MPI halo exchange
    io_restart_derivedtype.F90    # binary restart I/O
    fortran_utils.F90             # string utilities

    MOD_READ_BINARY_ARRAYS.F90    # binary array reader
    MOD_WRITE_BINARY_ARRAYS.F90   # binary array writer

    tracer_c_interface.F90    # C-compatible API (optional)
    tracer_init_from_mesh.F90 # tracer initialization from mesh files

  cmake/                  # CMake modules
  visualization/          # plotting scripts (plot_mesh.py, plot_scalars.py, plot_compare_precision.py)
  python/                 # Python scripts (C interface demos, test scripts)
  docs/                   # additional documentation
```
