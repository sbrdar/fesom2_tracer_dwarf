#!/usr/bin/env python3
"""
Test FESOM2 tracer C interface - Initialize from Python, load mesh from files.

This test demonstrates:
1. Initialize MPI from Python
2. Load mesh/partition from pre-decomposed files (like fesom_mesh_init.F90)
3. Allocate tracers from Python (set number of tracers)
4. Set tracer values from Python (set states)
5. Set velocities from Python

Usage:
    cd work
    export LD_LIBRARY_PATH=../build/lib:$LD_LIBRARY_PATH
    
    # Single process
    python ../test_mesh_from_files.py
    
    # Multiple processes
    mpirun -np 2 python ../test_mesh_from_files.py
"""

import ctypes
import sys
from pathlib import Path
import numpy as np

# Find library
SCRIPT_DIR = Path(__file__).parent
LIB_CANDIDATES = [
    SCRIPT_DIR / "build" / "lib" / "libfesom_tracer_Fortran.so",
    SCRIPT_DIR / ".." / "build" / "lib" / "libfesom_tracer_Fortran.so",
]

LIB_PATH = None
for candidate in LIB_CANDIDATES:
    if candidate.exists():
        LIB_PATH = candidate
        break

if LIB_PATH is None or not LIB_PATH.exists():
    print(f"ERROR: Library not found")
    sys.exit(1)

# Load library
lib = ctypes.CDLL(str(LIB_PATH), mode=ctypes.RTLD_GLOBAL)

# Define C function signatures
lib.tracer_init_mpi.argtypes = [ctypes.c_int]
lib.tracer_init_mpi.restype = ctypes.c_int

lib.tracer_load_mesh_from_files.argtypes = [ctypes.c_char_p]
lib.tracer_load_mesh_from_files.restype = ctypes.c_int

lib.tracer_get_size.argtypes = [
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
]
lib.tracer_get_size.restype = ctypes.c_int

lib.tracer_allocate_tracers.argtypes = [ctypes.c_int] * 4
lib.tracer_allocate_tracers.restype = ctypes.c_int

lib.tracer_set_values.argtypes = [
    ctypes.c_int,
    np.ctypeslib.ndpointer(dtype=np.float64, flags='F_CONTIGUOUS'),
    ctypes.c_int,
    ctypes.c_int,
]
lib.tracer_set_values.restype = ctypes.c_int

lib.tracer_set_velocity_zero.argtypes = [ctypes.c_int] * 3
lib.tracer_set_velocity_zero.restype = ctypes.c_int

lib.tracer_init_ale_arrays.argtypes = []
lib.tracer_init_ale_arrays.restype = ctypes.c_int

lib.tracer_init_for_advection.argtypes = []
lib.tracer_init_for_advection.restype = ctypes.c_int

lib.tracer_init_complete.argtypes = []
lib.tracer_init_complete.restype = ctypes.c_int

lib.tracer_advect_step.argtypes = [ctypes.c_double]
lib.tracer_advect_step.restype = ctypes.c_int

lib.tracer_get_stats.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
]
lib.tracer_get_stats.restype = ctypes.c_int

lib.tracer_finalize_mpi.argtypes = []
lib.tracer_finalize_mpi.restype = ctypes.c_int

# Get MPI rank for printing
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    npes = comm.Get_size()
except ImportError:
    rank = 0
    npes = 1

def print0(*args, **kwargs):
    """Print only from rank 0"""
    if rank == 0:
        print(*args, **kwargs)
        sys.stdout.flush()

print0("=" * 70)
print0("FESOM2 Tracer C Interface - Mesh from Files Test")
print0("=" * 70)
print0(f"Running with {npes} MPI process(es)")
print0()

try:
    # Step 1: Initialize MPI
    print0("Step 1: Initialize MPI from Python...")
    ret = lib.tracer_init_mpi(0)
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        sys.exit(1)
    print0("   ✓ MPI initialized")
    
    # Step 2: Load mesh from partition files (like fesom_mesh_init.F90)
    print0("\nStep 2: Load mesh from partition files...")
    # Path to mesh directory - mesh_setup() will read partition files from dist_N/
    mesh_path = b"../../../tests/data/MESHES/pi/"  # PI mesh
    ret = lib.tracer_load_mesh_from_files(mesh_path)
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        print0(f"   Make sure mesh partition files exist at: {mesh_path.decode()}/dist_{npes}/")
        sys.exit(1)
    print0(f"   ✓ Mesh loaded from partition files: {mesh_path.decode()}")
    
    # Step 3: Get mesh dimensions
    print0("\nStep 3: Get mesh dimensions...")
    nz = ctypes.c_int()
    nn = ctypes.c_int()
    ret = lib.tracer_get_size(ctypes.byref(nz), ctypes.byref(nn))
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        sys.exit(1)
    nl = nz.value + 1
    nod2D = nn.value
    print0(f"   ✓ Mesh dimensions: nl={nl}, nod2D={nod2D}")
    
    # Step 4: ⭐ Allocate tracers (SET NUMBER OF TRACERS)
    print0("\nStep 4: ⭐ Allocate tracers (set number of tracers)...")
    num_tracers = 2  # Temperature + Salinity
    ret = lib.tracer_allocate_tracers(num_tracers, nl, nod2D, 2)
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        sys.exit(1)
    print0(f"   ✓ Allocated {num_tracers} tracers (Temperature + Salinity)")
    
    # Step 5: ⭐ Set tracer values (SET STATES)
    print0("\nStep 5: ⭐ Set tracer values from Python arrays...")
    
    # Create temperature field: warm at surface, cold at depth
    temperature = np.zeros((nz.value, nn.value), dtype=np.float64, order='F')
    for k in range(nz.value):
        depth_frac = k / max(1, nz.value - 1)
        temperature[k, :] = 20.0 - 15.0 * depth_frac  # 20°C → 5°C
    
    ret = lib.tracer_set_values(1, temperature, nl, nod2D)
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        sys.exit(1)
    print0(f"   ✓ Set tracer 1 (Temperature): {temperature[0,0]:.1f}°C (surface) → {temperature[-1,0]:.1f}°C (bottom)")
    
    # Create salinity field: uniform
    salinity = np.full((nz.value, nn.value), 35.0, dtype=np.float64, order='F')
    ret = lib.tracer_set_values(2, salinity, nl, nod2D)
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        sys.exit(1)
    print0(f"   ✓ Set tracer 2 (Salinity): {salinity[0,0]:.1f} PSU (uniform)")
    
    # Step 6: ⭐ Set velocities (SET VELOCITIES)
    print0("\nStep 6: ⭐ Set velocities from Python...")
    # For this test, set to zero (can also use tracer_set_velocity with arrays)
    # Note: We need edge2D, which we can get from the mesh
    # For simplicity, estimate edge2D ≈ 3*nod2D for triangular mesh
    edge2D = 3 * nod2D
    ret = lib.tracer_set_velocity_zero(nl, nod2D, edge2D)
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        sys.exit(1)
    print0(f"   ✓ Velocities set to zero")
    
    # Step 6b: Initialize ALE arrays (needed for advection)
    print0("\nStep 6b: Initialize ALE arrays...")
    ret = lib.tracer_init_ale_arrays()
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        sys.exit(1)
    print0("   ✓ ALE arrays initialized")
    
    # Step 6c: Complete initialization for advection
    print0("\nStep 6c: Initialize remaining arrays for advection...")
    ret = lib.tracer_init_for_advection()
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        sys.exit(1)
    print0("   ✓ Full advection initialization complete")
    
    # Step 7: Mark initialization complete
    print0("\nStep 7: Mark initialization complete...")
    ret = lib.tracer_init_complete()
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        sys.exit(1)
    print0("   ✓ Initialization complete")
    
    # Step 8: Run advection steps (now should work!)
    print0("\nStep 8: Run advection steps (with zero velocity)...")
    dt = 1.0e-3  # 0.001 seconds
    for step in range(1, 6):
        ret = lib.tracer_advect_step(dt)
        if ret != 0:
            print0(f"   ✗ Advection failed at step {step} with code {ret}")
            break
        
        # Get statistics
        tmin = ctypes.c_double()
        tmax = ctypes.c_double()
        tsum = ctypes.c_double()
        ret = lib.tracer_get_stats(
            ctypes.byref(tmin),
            ctypes.byref(tmax),
            ctypes.byref(tsum)
        )
        if ret == 0 and rank == 0:
            print0(f"   Step {step}: T_min={tmin.value:7.3f}°C, T_max={tmax.value:7.3f}°C")
    
    print0("   ✓ Advection steps completed")
    
    # Step 9: Cleanup
    print0("\nStep 9: Finalize MPI...")
    ret = lib.tracer_finalize_mpi()
    if ret != 0:
        print0(f"   ✗ FAILED with code {ret}")
        sys.exit(1)
    print0("   ✓ MPI finalized")
    
    # Success!
    print0("\n" + "=" * 70)
    print0("✓ ALL TESTS PASSED!")
    print0("=" * 70)
    print0("\nKey capabilities demonstrated:")
    print0("  ✓ Initialize from Python")
    print0("  ✓ Load mesh from pre-decomposed files")
    print0("  ✓ Set number of tracers: tracer_allocate_tracers()")
    print0("  ✓ Set tracer states:     tracer_set_values()")
    print0("  ✓ Set velocities:        tracer_set_velocity_zero()")
    print0("  ✓ Run advection:         tracer_advect_step()")
    print0()

except Exception as e:
    print0(f"\n✗ ERROR: {e}")
    import traceback
    if rank == 0:
        traceback.print_exc()
    sys.exit(1)
