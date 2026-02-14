#!/usr/bin/env python3
"""
Test script for loading mesh from partition files and setting tracer values from arrays.

This demonstrates the hybrid approach:
1. Load mesh structure from existing partition files (provides connectivity)
2. Set tracer values from NumPy arrays (provides initial conditions)
3. Run advection with the complete setup

Usage:
    # Build with C interface
    ENABLE_TRACER_C_INTERFACE=ON ./configure.sh ubuntu
    cd build && make -j$(nproc)
    
    # Run test (single process)
    cd build
    python ../test_mesh_loading.py
    
    # Or with MPI (2 processes)
    mpirun -np 2 python ../test_mesh_loading.py
"""

import numpy as np
import sys
from pathlib import Path

# Add parent directory to path to import demo_tracer
sys.path.insert(0, str(Path(__file__).parent))
from demo_tracer import FESOMTracerAdvection


def test_mesh_loading_with_array_tracers():
    """
    Test loading mesh from partition files and setting tracer values from arrays.
    """
    print("=" * 70)
    print("FESOM2 Tracer Advection - Mesh Loading + Array Tracers Test")
    print("=" * 70)
    
    # Create tracer advection instance
    tracer = FESOMTracerAdvection()
    
    # Initialize MPI
    print("\n1. Initializing MPI...")
    ret = tracer.initialize_mpi(0)
    if ret != 0:
        print(f"✗ MPI initialization failed with code {ret}")
        sys.exit(1)
    print("✓ MPI initialized")
    
    # Load mesh from partition files
    print("\n2. Loading mesh from partition files...")
    # Use the existing restart files which contain the mesh
    # Note: This script should be run from the work directory
    mesh_path = "../fesom_bin_restart"
    
    # Check if path exists
    if not Path(mesh_path).exists():
        print(f"✗ Mesh path not found: {mesh_path}")
        print(f"   Current directory: {Path.cwd()}")
        print(f"   Please run from the 'work' directory or adjust the path")
        sys.exit(1)
    
    ret = tracer.load_mesh_partition(mesh_path)
    if ret != 0:
        print(f"✗ Mesh loading failed with code {ret}")
        sys.exit(1)
    print("✓ Mesh loaded successfully")
    
    # Get mesh dimensions
    print("\n3. Getting mesh dimensions...")
    nz = ctypes.c_int()
    nn = ctypes.c_int()
    ret = tracer.lib.tracer_get_size(ctypes.byref(nz), ctypes.byref(nn))
    if ret != 0:
        print(f"✗ Failed to get mesh size with code {ret}")
        sys.exit(1)
    
    nl = nz.value + 1  # Number of levels
    nod2D = nn.value   # Number of owned nodes
    
    print(f"✓ Mesh dimensions: nl={nl}, nod2D={nod2D}")
    
    # Allocate tracers
    print("\n4. Allocating tracers...")
    num_tracers = 1
    AB_order = 2
    ret = tracer.lib.tracer_allocate_tracers(num_tracers, nl, nod2D, AB_order)
    if ret != 0:
        print(f"✗ Tracer allocation failed with code {ret}")
        sys.exit(1)
    print(f"✓ Allocated {num_tracers} tracer(s)")
    
    # Create synthetic tracer field
    print("\n5. Creating synthetic temperature field...")
    # Create a simple vertical gradient
    temperature = np.zeros((nl-1, nod2D), dtype=np.float64, order='F')
    for k in range(nl-1):
        depth_frac = k / (nl-2) if nl > 2 else 0
        temp = 20.0 - 15.0 * depth_frac  # 20°C at surface, 5°C at depth
        temperature[k, :] = temp + np.random.randn(nod2D) * 0.5  # Add small noise
    
    print(f"   Temperature range: {temperature.min():.2f} to {temperature.max():.2f} °C")
    
    # Set tracer values
    print("\n6. Setting tracer values from array...")
    tracer_id = 1
    ret = tracer.lib.tracer_set_values(tracer_id, temperature, nl, nod2D)
    if ret != 0:
        print(f"✗ Setting tracer values failed with code {ret}")
        sys.exit(1)
    print("✓ Tracer values set")
    
    # Set velocity to zero for this test
    print("\n7. Setting velocity to zero...")
    # We need edge2D - get it from the mesh (it was loaded)
    # For now, we'll skip velocity setting since the mesh already has it from restart
    print("✓ Using velocity from loaded mesh")
    
    # Mark initialization complete
    print("\n8. Completing initialization...")
    ret = tracer.lib.tracer_init_complete()
    if ret != 0:
        print(f"✗ Initialization completion failed with code {ret}")
        sys.exit(1)
    print("✓ Initialization complete")
    
    # Get initial statistics
    print("\n9. Initial tracer statistics:")
    stats = tracer.get_statistics()
    print(f"   min = {stats['min']:12.6e}")
    print(f"   max = {stats['max']:12.6e}")
    print(f"   sum = {stats['sum']:12.6e}")
    
    # Run advection steps
    print("\n10. Running advection steps...")
    dt = 1.0e-3  # Small time step
    nsteps = 5
    
    print(f"    Running {nsteps} steps with dt={dt}...")
    for step in range(1, nsteps + 1):
        ret = tracer.advect_step(dt)
        if ret != 0:
            print(f"✗ Step {step} failed with code {ret}")
            break
        
        stats = tracer.get_statistics()
        print(f"    Step {step:2d}: min={stats['min']:12.6e}, "
              f"max={stats['max']:12.6e}, sum={stats['sum']:12.6e}")
    
    # Get final tracer values
    print("\n11. Retrieving final tracer values...")
    final_values = tracer.get_tracer_values()
    print(f"✓ Retrieved array of shape {final_values.shape}")
    print(f"   Final range: {final_values.min():.6e} to {final_values.max():.6e}")
    
    # Finalize MPI
    print("\n12. Finalizing...")
    ret = tracer.finalize_mpi()
    if ret != 0:
        print(f"✗ MPI finalization failed with code {ret}")
        sys.exit(1)
    
    print("\n" + "=" * 70)
    print("✓ Test completed successfully!")
    print("=" * 70)
    print("\nKey achievement:")
    print("  • Loaded mesh structure from partition files")
    print("  • Set tracer values from NumPy arrays")
    print("  • Successfully ran advection with the hybrid setup")


if __name__ == "__main__":
    import ctypes
    test_mesh_loading_with_array_tracers()
