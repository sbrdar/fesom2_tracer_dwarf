#!/usr/bin/env python3
"""
Simple MPI test for FESOM2 tracer C interface.
Tests the three key capabilities with minimal overhead.

Usage:
    cd work
    export LD_LIBRARY_PATH=../build/lib:$LD_LIBRARY_PATH
    mpirun -np 2 python ../test_mpi_simple.py
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

# Define signatures
lib.tracer_init_mpi.argtypes = [ctypes.c_int]
lib.tracer_init_mpi.restype = ctypes.c_int

lib.tracer_set_mesh_dims.argtypes = [ctypes.c_int] * 4
lib.tracer_set_mesh_dims.restype = ctypes.c_int

lib.tracer_set_partition.argtypes = [ctypes.c_int, ctypes.c_int]
lib.tracer_set_partition.restype = ctypes.c_int

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

lib.tracer_finalize_mpi.argtypes = []
lib.tracer_finalize_mpi.restype = ctypes.c_int

# Get MPI info
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    npes = comm.Get_size()
except ImportError:
    rank = 0
    npes = 1

# Only rank 0 prints
def print0(*args, **kwargs):
    if rank == 0:
        print(*args, **kwargs)
        sys.stdout.flush()

print0("=" * 60)
print0(f"FESOM2 Tracer C Interface - MPI Test (np={npes})")
print0("=" * 60)

try:
    # 1. Init MPI
    print0("\n1. Initialize MPI...")
    ret = lib.tracer_init_mpi(0)
    if ret != 0:
        print0(f"   FAILED: {ret}")
        sys.exit(1)
    print0("   ✓ OK")
    
    # 2. Set mesh
    print0("\n2. Set mesh dimensions...")
    nl, nod2D, edge2D, elem2D = 10, 5, 8, 4
    ret = lib.tracer_set_mesh_dims(nl, nod2D, edge2D, elem2D)
    if ret != 0:
        print0(f"   FAILED: {ret}")
        sys.exit(1)
    print0("   ✓ OK")
    
    # 3. Set partition (SKIP - needs mesh connectivity)
    print0("\n3. Set partition... SKIPPED (needs mesh connectivity)")
    # Note: tracer_set_partition calls init_mpi_types which needs mesh arrays
    # For this simple test, we skip it and allocate with nod2D directly
    
    # 4. Allocate tracers (KEY!)
    print0("\n4. ⭐ Allocate 2 tracers...")
    ret = lib.tracer_allocate_tracers(2, nl, nod2D, 2)
    if ret != 0:
        print0(f"   FAILED: {ret}")
        sys.exit(1)
    print0("   ✓ OK - 2 tracers allocated")
    
    # 5. Set tracer values (KEY!)
    print0("\n5. ⭐ Set tracer values...")
    temp = np.full((nl-1, nod2D), 15.0, dtype=np.float64, order='F')
    ret = lib.tracer_set_values(1, temp, nl, nod2D)
    if ret != 0:
        print0(f"   FAILED: {ret}")
        sys.exit(1)
    print0("   ✓ OK - Temperature set to 15°C")
    
    salt = np.full((nl-1, nod2D), 35.0, dtype=np.float64, order='F')
    ret = lib.tracer_set_values(2, salt, nl, nod2D)
    if ret != 0:
        print0(f"   FAILED: {ret}")
        sys.exit(1)
    print0("   ✓ OK - Salinity set to 35 PSU")
    
    # 6. Set velocity (KEY!)
    print0("\n6. ⭐ Set velocity to zero...")
    ret = lib.tracer_set_velocity_zero(nl, nod2D, edge2D)
    if ret != 0:
        print0(f"   FAILED: {ret}")
        sys.exit(1)
    print0("   ✓ OK - Velocity = 0")
    
    # 7. Cleanup
    print0("\n7. Finalize MPI...")
    ret = lib.tracer_finalize_mpi()
    if ret != 0:
        print0(f"   FAILED: {ret}")
        sys.exit(1)
    print0("   ✓ OK")
    
    print0("\n" + "=" * 60)
    print0("✓ ALL TESTS PASSED!")
    print0("=" * 60)
    print0("\nKey capabilities verified:")
    print0("  ✓ Set number of tracers: tracer_allocate_tracers()")
    print0("  ✓ Set tracer states:     tracer_set_values()")
    print0("  ✓ Set velocities:        tracer_set_velocity_zero()")
    print0()

except Exception as e:
    print0(f"\nERROR: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
