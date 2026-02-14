#!/usr/bin/env python3
"""Load and plot FESOM2 mesh from raw binary files written by --save-mesh.

Usage:
    python plot_mesh.py [mesh_output_dir]

Default mesh_output_dir is 'mesh_output'.
Saves mesh_plot.png in the mesh directory and shows the plot interactively.
"""

import sys
import os
import numpy as np


def load_mesh(dirpath):
    """Read mesh_info.txt and binary files, return arrays."""
    info = {}
    with open(os.path.join(dirpath, "mesh_info.txt")) as f:
        for line in f:
            key, val = line.strip().split("=")
            info[key] = int(val)

    nod2D = info["nod2D"]
    elem2D = info["elem2D"]
    nl = info["nl"]

    # Fortran (2, nod2D) column-major -> numpy (2, nod2D) with order='F'
    coords = np.fromfile(
        os.path.join(dirpath, "coord_nod2D.bin"), dtype=np.float64
    ).reshape((2, nod2D), order="F")
    x = coords[0, :]
    y = coords[1, :]

    # Fortran (3, elem2D) column-major, 0-based indices
    triangles = np.fromfile(
        os.path.join(dirpath, "elem2D_nodes.bin"), dtype=np.int32
    ).reshape((3, elem2D), order="F").T  # matplotlib wants (elem2D, 3)

    zbar = np.fromfile(os.path.join(dirpath, "zbar.bin"), dtype=np.float64)
    depth = np.fromfile(os.path.join(dirpath, "depth.bin"), dtype=np.float64)

    return x, y, triangles, zbar, depth


def plot_surface_mesh(x, y, triangles, save_path=None):
    """Plot the surface triangulation."""
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri

    triang = tri.Triangulation(x, y, triangles)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8))
    ax.triplot(triang, linewidth=0.5, color="steelblue")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    ax.set_title(f"Surface mesh: {len(x)} nodes, {len(triangles)} triangles")
    ax.set_aspect("equal")

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"Saved plot to {save_path}")

    plt.show()


if __name__ == "__main__":
    mesh_dir = sys.argv[1] if len(sys.argv) > 1 else "mesh_output"

    if not os.path.isdir(mesh_dir):
        print(f"Error: directory '{mesh_dir}' not found", file=sys.stderr)
        sys.exit(1)

    x, y, triangles, zbar, depth = load_mesh(mesh_dir)

    print(f"Loaded mesh from '{mesh_dir}':")
    print(f"  Nodes:    {len(x)}")
    print(f"  Elements: {len(triangles)}")
    print(f"  Levels:   {len(zbar)}")
    print(f"  x range:  [{x.min():.1f}, {x.max():.1f}] m")
    print(f"  y range:  [{y.min():.1f}, {y.max():.1f}] m")
    print(f"  depth:    [{depth.min():.1f}, {depth.max():.1f}] m")
    print(f"  zbar:     [{zbar.min():.1f}, {zbar.max():.1f}] m")

    plot_surface_mesh(x, y, triangles, save_path=os.path.join(mesh_dir, "mesh_plot.png"))
