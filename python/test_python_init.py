#!/usr/bin/env python3
"""
Test script for Python array initialization of FESOM2 tracer advection.

This demonstrates initializing the tracer advection from NumPy arrays
without requiring binary restart files.

Usage:
    # Build with C interface
    ENABLE_TRACER_C_INTERFACE=ON ./configure.sh ubuntu
    cd build && make -j$(nproc)
    
    # Run test
    cd work
    python ../test_python_init.py
"""

import numpy as np
import sys
from pathlib import Path

# Add parent directory to path to import demo_tracer
sys.path.insert(0, str(Path(__file__).parent))
from demo_tracer import FESOMTracerAdvection


def create_simple_test_case():
    """
    Create a simple synthetic test case for tracer advection.
    
    Returns:
        dict with mesh_dims, levels, thickness, partition, tracer_values
    """
    # Simple 1D column test: 10 vertical levels, 5 nodes
    nl = 10
    nod2D = 5
    edge2D = 8  # Approximate for simple mesh
    elem2D = 4
    
    # All nodes have same vertical structure
    nlevels = np.full(nod2D, nl-1, dtype=np.int32)  # Bottom at level 9
    ulevels = np.ones(nod2D, dtype=np.int32)        # Surface at level 1
    
    # Uniform layer thickness: 100m per layer
    thickness = np.ones((nl-1, nod2D), dtype=np.float64, order='F') * 100.0
    
    # Create temperature field with vertical gradient
    # Warm at surface (20°C), cold at depth (5°C)
    temperature = np.zeros((nl-1, nod2D), dtype=np.float64, order='F')
    for k in range(nl-1):
        depth_frac = k / (nl-2)  # 0 at surface, 1 at bottom
        temp = 20.0 - 15.0 * depth_frac  # Linear gradient
        temperature[k, :] = temp + np.random.randn(nod2D) * 0.1  # Add small noise
    
    return {
        'mesh_dims': {
            'nl': nl,
            'nod2D': nod2D,
            'edge2D': edge2D,
            'elem2D': elem2D
        },
        'levels': {
            'nlevels': nlevels,
            'ulevels': ulevels
        },
        'thickness': thickness,
        'partition': {
            'myDim_nod2D': nod2D,
            'eDim_nod2D': 0  # No halo for single process
        },
        'tracer_values': {
            1: temperature  # Tracer ID 1 = temperature
        }
    }


def test_array_initialization():
    """
    Test initialization from NumPy arrays.
    """
    print("=" * 70)
    print("FESOM2 Tracer Advection - Array Initialization Test")
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
    
    # Create test case
    print("\n2. Creating synthetic test case...")
    test_case = create_simple_test_case()
    print(f"✓ Created test case:")
    print(f"   - Mesh: {test_case['mesh_dims']['nl']} levels, "
          f"{test_case['mesh_dims']['nod2D']} nodes")
    print(f"   - Temperature range: "
          f"{test_case['tracer_values'][1].min():.2f} to "
          f"{test_case['tracer_values'][1].max():.2f} °C")
    
    # Initialize from arrays
    print("\n3. Initializing from NumPy arrays...")
    ret = tracer.initialize_from_arrays(
        mesh_dims=test_case['mesh_dims'],
        levels=test_case['levels'],
        thickness=test_case['thickness'],
        partition=test_case['partition'],
        tracer_values=test_case['tracer_values'],
        velocity=None  # Zero velocity for this test
    )
    if ret != 0:
        print(f"✗ Array initialization failed with code {ret}")
        sys.exit(1)
    
    # Get initial statistics
    print("\n4. Initial tracer statistics:")
    stats = tracer.get_statistics()
    print(f"   min = {stats['min']:12.6e}")
    print(f"   max = {stats['max']:12.6e}")
    print(f"   sum = {stats['sum']:12.6e}")
    
    # Run advection steps (with zero velocity, should be stable)
    print("\n5. Running advection steps (zero velocity)...")
    dt = 1.0e-3  # Small time step
    nsteps = 5
    
    for step in range(1, nsteps + 1):
        ret = tracer.advect_step(dt)
        if ret != 0:
            print(f"✗ Step {step} failed with code {ret}")
            break
        
        stats = tracer.get_statistics()
        print(f"   Step {step:2d}: min={stats['min']:12.6e}, "
              f"max={stats['max']:12.6e}, sum={stats['sum']:12.6e}")
    
    # Get final tracer values
    print("\n6. Retrieving final tracer values...")
    final_values = tracer.get_tracer_values()
    print(f"✓ Retrieved array of shape {final_values.shape}")
    print(f"   Final range: {final_values.min():.6e} to {final_values.max():.6e}")
    
    # Check conservation (with zero velocity, tracer should be conserved)
    initial_sum = test_case['tracer_values'][1].sum()
    final_sum = final_values.sum()
    rel_change = abs(final_sum - initial_sum) / abs(initial_sum)
    print(f"\n7. Conservation check:")
    print(f"   Initial sum: {initial_sum:.6e}")
    print(f"   Final sum:   {final_sum:.6e}")
    print(f"   Relative change: {rel_change:.6e}")
    if rel_change < 1e-10:
        print("   ✓ Excellent conservation!")
    elif rel_change < 1e-6:
        print("   ✓ Good conservation")
    else:
        print("   ⚠ Warning: significant change in tracer sum")
    
    # Finalize MPI
    print("\n8. Finalizing...")
    ret = tracer.finalize_mpi()
    if ret != 0:
        print(f"✗ MPI finalization failed with code {ret}")
        sys.exit(1)
    
    print("\n" + "=" * 70)
    print("✓ Test completed successfully!")
    print("=" * 70)


def test_with_velocity():
    """
    Test with non-zero velocity field (future enhancement).
    """
    print("\n" + "=" * 70)
    print("Test with velocity field - NOT YET IMPLEMENTED")
    print("=" * 70)
    print("To implement this test:")
    print("1. Create realistic velocity field (uv, w)")
    print("2. Pass velocity dict to initialize_from_arrays()")
    print("3. Observe tracer advection")


if __name__ == "__main__":
    # Run basic test
    test_array_initialization()
    
    # Future: test with velocity
    # test_with_velocity()
