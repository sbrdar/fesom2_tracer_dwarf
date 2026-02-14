#!/usr/bin/env python3
"""
Python Demo for FESOM2 Tracer Advection Library

This script demonstrates how to use the libfesom_tracer_Fortran.so shared library
from Python using ctypes. It replicates the functionality of the fesom.F90 program.

NOTE: This requires building with ENABLE_TRACER_C_INTERFACE=ON and implementing
the tracer_c_interface.F90 wrapper module.

Usage:
    # Build with C interface
    ENABLE_TRACER_C_INTERFACE=ON ./configure.sh ubuntu
    cd build && make -j$(nproc)
    
    # Run demo
    cd work
    python ../demo_tracer.py
"""

import ctypes
import os
import sys
from pathlib import Path
import numpy as np

# Path to the built library
BUILD_DIR = Path(__file__).parent / "lib"
TRACER_LIB = BUILD_DIR / "libfesom_tracer_Fortran.so"

# Check if library exists
if not TRACER_LIB.exists():
    print(f"Error: Tracer library not found at {TRACER_LIB}")
    print("Please build the project first with C interface enabled:")
    print("  ENABLE_TRACER_C_INTERFACE=ON ./configure.sh ubuntu")
    print("  cd build && make -j$(nproc)")
    sys.exit(1)


class FESOMTracerAdvection:
    """
    Python wrapper for FESOM tracer advection library.
    
    This class provides a high-level interface to the tracer advection
    functionality implemented in Fortran.
    
    NOTE: Requires C interface implementation in tracer_c_interface.F90
    """
    
    def __init__(self):
        """Initialize the library handle."""
        # Load library
        self.lib = ctypes.CDLL(str(TRACER_LIB), mode=ctypes.RTLD_GLOBAL)
        print(f"✓ Loaded tracer library: {TRACER_LIB}")
        
        # Define C function signatures
        
        # tracer_init_mpi
        self.lib.tracer_init_mpi.argtypes = [ctypes.c_int]  # mpi_comm
        self.lib.tracer_init_mpi.restype = ctypes.c_int
        
        # tracer_init
        self.lib.tracer_init.argtypes = [ctypes.c_char_p]  # restart_dir
        self.lib.tracer_init.restype = ctypes.c_int
        
        # tracer_advect_step
        self.lib.tracer_advect_step.argtypes = [ctypes.c_double]  # dt
        self.lib.tracer_advect_step.restype = ctypes.c_int
        
        # tracer_get_stats
        self.lib.tracer_get_stats.argtypes = [
            ctypes.POINTER(ctypes.c_double),  # tmin
            ctypes.POINTER(ctypes.c_double),  # tmax
            ctypes.POINTER(ctypes.c_double),  # tsum
        ]
        self.lib.tracer_get_stats.restype = ctypes.c_int
        
        # tracer_get_size
        self.lib.tracer_get_size.argtypes = [
            ctypes.POINTER(ctypes.c_int),  # nz
            ctypes.POINTER(ctypes.c_int),  # nn
        ]
        self.lib.tracer_get_size.restype = ctypes.c_int
        
        # tracer_get_values
        self.lib.tracer_get_values.argtypes = [
            np.ctypeslib.ndpointer(dtype=np.float64, flags='F_CONTIGUOUS'),  # values
            ctypes.c_int,  # nz
            ctypes.c_int,  # nn
        ]
        self.lib.tracer_get_values.restype = ctypes.c_int
        
        # tracer_finalize
        self.lib.tracer_finalize.argtypes = [ctypes.c_char_p]  # output_dir
        self.lib.tracer_finalize.restype = ctypes.c_int
        
        # tracer_finalize_mpi
        self.lib.tracer_finalize_mpi.argtypes = []
        self.lib.tracer_finalize_mpi.restype = ctypes.c_int
        
    def initialize_mpi(self, mpi_comm=0):
        """
        Initialize MPI (if not already initialized).
        
        Args:
            mpi_comm: MPI communicator (0 = MPI_COMM_WORLD)
        
        Returns:
            0 on success, non-zero on error
        """
        return self.lib.tracer_init_mpi(mpi_comm)
    
    def initialize(self, restart_path="./fesom_bin_restart"):
        """
        Initialize tracer advection from binary restart files.
        
        Args:
            restart_path: Path to restart directory
        
        Returns:
            0 on success, non-zero on error
        """
        print(f"Initializing from: {restart_path}")
        return self.lib.tracer_init(restart_path.encode())
    
    def advect_step(self, dt=1.0e-3):
        """
        Perform one tracer advection time step.
        
        Args:
            dt: Time step size
        
        Returns:
            0 on success, non-zero on error
        """
        return self.lib.tracer_advect_step(dt)
    
    def get_tracer_values(self):
        """
        Get current tracer values as numpy array.
        
        Returns:
            numpy array of tracer values (nz x nn)
        """
        # Get array size
        nz = ctypes.c_int()
        nn = ctypes.c_int()
        ret = self.lib.tracer_get_size(ctypes.byref(nz), ctypes.byref(nn))
        if ret != 0:
            raise RuntimeError(f"tracer_get_size failed with code {ret}")
        
        # Allocate array (Fortran column-major order)
        values = np.zeros((nz.value, nn.value), dtype=np.float64, order='F')
        
        # Get values
        ret = self.lib.tracer_get_values(values, nz, nn)
        if ret != 0:
            raise RuntimeError(f"tracer_get_values failed with code {ret}")
        
        return values
    
    def get_statistics(self):
        """
        Get min/max/sum statistics of tracer values.
        
        Returns:
            dict with 'min', 'max', 'sum' keys
        """
        tmin = ctypes.c_double()
        tmax = ctypes.c_double()
        tsum = ctypes.c_double()
        
        ret = self.lib.tracer_get_stats(
            ctypes.byref(tmin),
            ctypes.byref(tmax),
            ctypes.byref(tsum)
        )
        
        if ret != 0:
            raise RuntimeError(f"tracer_get_stats failed with code {ret}")
        
        return {
            'min': tmin.value,
            'max': tmax.value,
            'sum': tsum.value
        }
    
    def finalize(self, output_path="./fesom_bin_restart"):
        """
        Write updated restart files and cleanup.
        
        Args:
            output_path: Path to write restart files
        
        Returns:
            0 on success, non-zero on error
        """
        print(f"Writing restart files to: {output_path}")
        return self.lib.tracer_finalize(output_path.encode())
    
    def finalize_mpi(self):
        """
        Finalize MPI (if we initialized it).
        
        Returns:
            0 on success, non-zero on error
        """
        return self.lib.tracer_finalize_mpi()


def main():
    """
    Main demo function - replicates fesom.F90 functionality.
    """
    print("=" * 60)
    print("FESOM2 Tracer Advection - Python Demo")
    print("=" * 60)
    
    # Initialize tracer advection
    tracer = FESOMTracerAdvection()
    
    # Initialize MPI
    ret = tracer.initialize_mpi(0)  # 0 = MPI_COMM_WORLD
    if ret != 0:
        print(f"Error: MPI initialization failed with code {ret}")
        sys.exit(1)
    
    # Initialize from restart files
    restart_path = "./fesom_bin_restart"
    
    ret = tracer.initialize(restart_path)
    if ret != 0:
        print(f"Error: Initialization failed with code {ret}")
        sys.exit(1)
    
    print("\nRunning 10 advection steps...")
    
    # Run 10 time steps (like fesom.F90)
    dt = 1.0e-3
    for step in range(1, 11):
        ret = tracer.advect_step(dt)
        if ret != 0:
            print(f"Error: Step {step} failed with code {ret}")
            break
        
        # Get statistics
        stats = tracer.get_statistics()
        print(f"Step {step:2d}: min={stats['min']:12.6e}, "
              f"max={stats['max']:12.6e}, sum={stats['sum']:12.6e}")
    
    # Write output and cleanup
    ret = tracer.finalize(restart_path)
    if ret != 0:
        print(f"Error: Finalization failed with code {ret}")
        sys.exit(1)
    
    # Finalize MPI
    ret = tracer.finalize_mpi()
    if ret != 0:
        print(f"Error: MPI finalization failed with code {ret}")
        sys.exit(1)
    
    print("\n" + "=" * 60)
    print("Demo completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
