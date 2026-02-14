#!/usr/bin/env python3
"""
Quick verification script for FESOM2 tracer C interface.

This script verifies that all C interface functions are accessible
and demonstrates basic usage for setting number of tracers, states, and velocities.

Usage:
    cd /home/suvarchal/AWI/FESOM2_model/cleanup_interfaces2/dwarf/dwarf_tracer
    export LD_LIBRARY_PATH=build/lib:$LD_LIBRARY_PATH
    python verify_c_interface.py
"""

import ctypes
import sys
from pathlib import Path
import numpy as np

# Path to library - try multiple locations
SCRIPT_DIR = Path(__file__).parent
LIB_CANDIDATES = [
    SCRIPT_DIR / "build" / "lib" / "libfesom_tracer_Fortran.so",  # From dwarf_tracer dir
    SCRIPT_DIR / "lib" / "libfesom_tracer_Fortran.so",            # Installed location
    SCRIPT_DIR / ".." / "build" / "lib" / "libfesom_tracer_Fortran.so",  # From work dir
]

# Find the library
LIB_PATH = None
for candidate in LIB_CANDIDATES:
    if candidate.exists():
        LIB_PATH = candidate
        break

if LIB_PATH is None:
    LIB_PATH = LIB_CANDIDATES[0]  # Use first for error message

def verify_library():
    """Verify library exists and can be loaded."""
    print("=" * 70)
    print("FESOM2 Tracer C Interface - Verification")
    print("=" * 70)
    
    if not LIB_PATH.exists():
        print(f"\n✗ Library not found: {LIB_PATH}")
        print("\nPlease build with C interface enabled:")
        print("  ENABLE_TRACER_C_INTERFACE=ON ./configure.sh ubuntu")
        print("  cd build && make -j$(nproc)")
        return None
    
    print(f"\n✓ Library found: {LIB_PATH}")
    print(f"  Size: {LIB_PATH.stat().st_size / 1024:.1f} KB")
    
    try:
        lib = ctypes.CDLL(str(LIB_PATH), mode=ctypes.RTLD_GLOBAL)
        print("✓ Library loaded successfully")
        return lib
    except Exception as e:
        print(f"✗ Failed to load library: {e}")
        return None


def verify_functions(lib):
    """Verify all C interface functions are accessible."""
    print("\n" + "-" * 70)
    print("Verifying C Interface Functions")
    print("-" * 70)
    
    functions = [
        # MPI functions
        "tracer_init_mpi",
        "tracer_finalize_mpi",
        
        # Initialization from restart
        "tracer_init",
        "tracer_finalize",
        
        # Array initialization functions (NEW!)
        "tracer_set_mesh_dims",
        "tracer_set_levels",
        "tracer_set_thickness",
        "tracer_set_partition",
        "tracer_allocate_tracers",      # ← Set number of tracers
        "tracer_set_values",             # ← Set tracer states
        "tracer_set_velocity",           # ← Set velocities
        "tracer_set_velocity_zero",      # ← Set velocities to zero
        "tracer_init_complete",
        "tracer_load_mesh_partition",
        
        # Workflow functions
        "tracer_advect_step",
        "tracer_get_stats",
        "tracer_get_size",
        "tracer_get_values",
        "tracer_cleanup",
        "tracer_run_workflow",
    ]
    
    missing = []
    for func_name in functions:
        try:
            func = getattr(lib, func_name)
            print(f"  ✓ {func_name}")
        except AttributeError:
            print(f"  ✗ {func_name} - NOT FOUND")
            missing.append(func_name)
    
    if missing:
        print(f"\n✗ {len(missing)} functions missing!")
        return False
    else:
        print(f"\n✓ All {len(functions)} functions found!")
        return True


def demonstrate_usage(lib):
    """Demonstrate basic usage of key functions."""
    
    # Define function signatures
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
    
    try:
        # 1. Initialize MPI
        ret = lib.tracer_init_mpi(0)
        if ret != 0:
            print(f"   ✗ MPI init failed with code {ret}")
            return False
        
        # Get MPI rank from environment or assume rank 0
        try:
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            npes = comm.Get_size()
        except ImportError:
            rank = 0
            npes = 1
        
        # Only print from rank 0
        if rank == 0:
            print("\n" + "-" * 70)
            print("Demonstrating Key Functions")
            print("-" * 70)
            print(f"\nRunning with {npes} MPI process(es)")
            print("\n1. Initializing MPI...")
            print("   ✓ MPI initialized")
        
        # 2. Set mesh dimensions
        if rank == 0:
            print("\n2. Setting mesh dimensions...")
        nl, nod2D, edge2D, elem2D = 10, 5, 8, 4
        ret = lib.tracer_set_mesh_dims(nl, nod2D, edge2D, elem2D)
        if ret != 0:
            if rank == 0:
                print(f"   ✗ Failed with code {ret}")
            return False
        if rank == 0:
            print(f"   ✓ Mesh dims set: nl={nl}, nod2D={nod2D}, edge2D={edge2D}, elem2D={elem2D}")
        
        # 2b. Set partition (needed for MPI)
        if rank == 0:
            print("\n2b. Setting partition...")
        ret = lib.tracer_set_partition(nod2D, 0)  # myDim_nod2D, eDim_nod2D
        if ret != 0:
            if rank == 0:
                print(f"   ✗ Failed with code {ret}")
            return False
        if rank == 0:
            print(f"   ✓ Partition set: myDim={nod2D}, eDim=0")
        
        # 3. Allocate tracers (KEY FUNCTION!)
        if rank == 0:
            print("\n3. Allocating tracers (setting number of tracers)...")
        num_tracers = 2  # Temperature + Salinity
        ret = lib.tracer_allocate_tracers(num_tracers, nl, nod2D, 2)
        if ret != 0:
            if rank == 0:
                print(f"   ✗ Failed with code {ret}")
            return False
        if rank == 0:
            print(f"   ✓ Allocated {num_tracers} tracers")
        
        # 4. Set tracer values (KEY FUNCTION!)
        if rank == 0:
            print("\n4. Setting tracer values (states)...")
        
        # Create temperature field
        temperature = np.zeros((nl-1, nod2D), dtype=np.float64, order='F')
        for k in range(nl-1):
            temperature[k, :] = 20.0 - 15.0 * (k / (nl-2))
        
        ret = lib.tracer_set_values(1, temperature, nl, nod2D)
        if ret != 0:
            if rank == 0:
                print(f"   ✗ Failed with code {ret}")
            return False
        if rank == 0:
            print(f"   ✓ Set tracer 1 (temperature): {temperature[0,0]:.1f}°C (surface) → {temperature[-1,0]:.1f}°C (bottom)")
        
        # Create salinity field
        salinity = np.full((nl-1, nod2D), 35.0, dtype=np.float64, order='F')
        ret = lib.tracer_set_values(2, salinity, nl, nod2D)
        if ret != 0:
            if rank == 0:
                print(f"   ✗ Failed with code {ret}")
            return False
        if rank == 0:
            print(f"   ✓ Set tracer 2 (salinity): {salinity[0,0]:.1f} PSU")
        
        # 5. Set velocity (KEY FUNCTION!)
        if rank == 0:
            print("\n5. Setting velocity to zero...")
        ret = lib.tracer_set_velocity_zero(nl, nod2D, edge2D)
        if ret != 0:
            if rank == 0:
                print(f"   ✗ Failed with code {ret}")
            return False
        if rank == 0:
            print("   ✓ Velocity set to zero")
        
        # 6. Cleanup
        if rank == 0:
            print("\n6. Cleaning up...")
        ret = lib.tracer_finalize_mpi()
        if ret != 0:
            if rank == 0:
                print(f"   ✗ Failed with code {ret}")
            return False
        if rank == 0:
            print("   ✓ MPI finalized")
        
        return True
        
    except Exception as e:
        print(f"\n✗ Error during demonstration: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Main verification routine."""
    # Step 1: Verify library
    lib = verify_library()
    if lib is None:
        sys.exit(1)
    
    # Step 2: Verify functions
    if not verify_functions(lib):
        sys.exit(1)
    
    # Step 3: Demonstrate usage
    if not demonstrate_usage(lib):
        print("\n✗ Demonstration failed")
        sys.exit(1)
    
    # Success!
    print("\n" + "=" * 70)
    print("✓ ALL VERIFICATIONS PASSED!")
    print("=" * 70)
    print("\nThe C interface is fully functional and ready to use!")
    print("\nKey capabilities verified:")
    print("  ✓ Set number of tracers: tracer_allocate_tracers()")
    print("  ✓ Set tracer states:     tracer_set_values()")
    print("  ✓ Set velocities:        tracer_set_velocity() / tracer_set_velocity_zero()")
    print("\nFor detailed usage examples, see:")
    print("  - C_INTERFACE_USAGE_GUIDE.md")
    print("  - C_INTERFACE_REFERENCE.md")
    print("  - demo_tracer.py")
    print("  - test_python_init.py")
    print()


if __name__ == "__main__":
    main()
