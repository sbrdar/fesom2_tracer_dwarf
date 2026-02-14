#!/usr/bin/env python3
"""Load and plot FESOM2 tracer scalar fields from raw binary files written by --save-scalars.

Usage:
    python plot_scalars.py [dir]                     # single-directory plot
    python plot_scalars.py --compare dir1 dir2       # side-by-side + difference

Default dir is 'mesh_output'.
Requires --save-mesh and --save-scalars to have been run together so that both
mesh geometry and scalar data are present in each output directory.

Saves per-step PNG images for each tracer's surface layer (k=0) and an animation GIF.
"""

import sys
import os
import argparse
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

    coords = np.fromfile(
        os.path.join(dirpath, "coord_nod2D.bin"), dtype=np.float64
    ).reshape((2, nod2D), order="F")
    x = coords[0, :]
    y = coords[1, :]

    triangles = np.fromfile(
        os.path.join(dirpath, "elem2D_nodes.bin"), dtype=np.int32
    ).reshape((3, elem2D), order="F").T

    zbar = np.fromfile(os.path.join(dirpath, "zbar.bin"), dtype=np.float64)
    depth = np.fromfile(os.path.join(dirpath, "depth.bin"), dtype=np.float64)

    return x, y, triangles, zbar, depth


def load_scalars(dirpath):
    """Read scalar_info.txt and binary tracer files, return info dict and tracer arrays.

    Returns:
        info: dict with keys nod2D, nz, nsteps, num_tracers, tracer_names, wp
        data: dict mapping tracer name -> np.ndarray of shape (nsteps+1, nod2D, nz)
    """
    info = {}
    tracer_names = []
    with open(os.path.join(dirpath, "scalar_info.txt")) as f:
        for line in f:
            key, val = line.strip().split("=", 1)
            if key.startswith("tracer_"):
                tracer_names.append(val)
            else:
                info[key] = int(val)

    nod2D = info["nod2D"]
    nz = info["nz"]
    nsteps = info["nsteps"]
    num_tracers = info["num_tracers"]
    nframes = nsteps + 1

    info["tracer_names"] = tracer_names

    data = {}
    for name in tracer_names:
        fpath = os.path.join(dirpath, f"{name}.bin")
        raw = np.fromfile(fpath, dtype=np.float64)
        expected = nframes * nz * nod2D
        if raw.size != expected:
            raise ValueError(
                f"{name}.bin: expected {expected} values "
                f"({nframes} frames x {nz} layers x {nod2D} nodes), got {raw.size}"
            )
        # Fortran column-major (nz, nod2D) per frame -> C (nod2D, nz)
        frames = raw.reshape((nframes, nod2D, nz))
        data[name] = frames  # shape: (nframes, nod2D, nz)

    return info, data


def make_triangulation(x, y, triangles):
    """Create Triangulation with wrap-around masking for periodic meshes."""
    import matplotlib.tri as tri

    triang = tri.Triangulation(x, y, triangles)

    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    if x_range > 0 and y_range > 0:
        tri_x = x[triangles]
        tri_y = y[triangles]
        dx_max = tri_x.max(axis=1) - tri_x.min(axis=1)
        dy_max = tri_y.max(axis=1) - tri_y.min(axis=1)
        wrap_mask = (dx_max > 0.5 * x_range) | (dy_max > 0.5 * y_range)
        if wrap_mask.any():
            triang.set_mask(wrap_mask)
            print(f"  Masked {wrap_mask.sum()} wrap-around triangles (periodic mesh)")

    return triang


def get_precision_label(info):
    """Return a short label like 'DP', 'SP', or 'HP' based on wp in info."""
    wp = info.get("wp", 8)
    if wp == 2:
        return "HP"
    elif wp == 4:
        return "SP"
    elif wp == 8:
        return "DP"
    else:
        return f"WP{wp}"


def plot_scalars(dirpath):
    """Plot surface layer of each tracer at each timestep."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x, y, triangles, zbar, depth = load_mesh(dirpath)
    info, data = load_scalars(dirpath)

    nsteps = info["nsteps"]
    nframes = nsteps + 1
    tracer_names = info["tracer_names"]
    prec_label = get_precision_label(info)

    triang = make_triangulation(x, y, triangles)

    img_dir = os.path.join(dirpath, "scalar_plots")
    os.makedirs(img_dir, exist_ok=True)

    for name in tracer_names:
        arr = data[name]  # (nframes, nod2D, nz)
        surface = arr[:, :, 0]

        vmin = surface.min()
        vmax = surface.max()
        if vmin == vmax:
            vmin -= 0.5
            vmax += 0.5

        frames_for_gif = []
        for step in range(nframes):
            fig, ax = plt.subplots(1, 1, figsize=(8, 7))
            tpc = ax.tripcolor(triang, surface[step], shading="flat",
                               cmap="RdYlBu_r", vmin=vmin, vmax=vmax)
            fig.colorbar(tpc, ax=ax, label=name.capitalize())
            ax.set_xlabel("x (m)")
            ax.set_ylabel("y (m)")
            ax.set_title(f"{name.capitalize()} [{prec_label}] surface (k=0) — step {step}/{nsteps}")
            ax.set_aspect("equal")

            fname = os.path.join(img_dir, f"{name}_step{step:04d}.png")
            fig.savefig(fname, dpi=100, bbox_inches="tight")
            frames_for_gif.append(fname)
            plt.close(fig)

        print(f"  Saved {nframes} frames for '{name}' to {img_dir}/")

        try:
            from PIL import Image
            imgs = [Image.open(f) for f in frames_for_gif]
            gif_path = os.path.join(dirpath, f"{name}_animation.gif")
            imgs[0].save(gif_path, save_all=True, append_images=imgs[1:],
                         duration=300, loop=0)
            print(f"  Animation saved to {gif_path}")
        except ImportError:
            print("  (Pillow not installed — skipping GIF animation)")


def plot_compare(dir1, dir2, output_dir=None):
    """Plot two runs side by side with difference panel."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x1, y1, tri1, _, _ = load_mesh(dir1)
    x2, y2, tri2, _, _ = load_mesh(dir2)
    info1, data1 = load_scalars(dir1)
    info2, data2 = load_scalars(dir2)

    label1 = get_precision_label(info1)
    label2 = get_precision_label(info2)

    # Validate compatible grids
    if info1["nod2D"] != info2["nod2D"]:
        raise ValueError(f"Node count mismatch: {info1['nod2D']} vs {info2['nod2D']}")
    if info1["nsteps"] != info2["nsteps"]:
        raise ValueError(f"Step count mismatch: {info1['nsteps']} vs {info2['nsteps']}")

    nsteps = info1["nsteps"]
    nframes = nsteps + 1
    tracer_names = info1["tracer_names"]

    # Use first directory's mesh for triangulation (same grid)
    triang = make_triangulation(x1, y1, tri1)

    if output_dir is None:
        output_dir = f"compare_{label1}_vs_{label2}"
    os.makedirs(output_dir, exist_ok=True)
    img_dir = os.path.join(output_dir, "compare_plots")
    os.makedirs(img_dir, exist_ok=True)

    for name in tracer_names:
        surf1 = data1[name][:, :, 0]  # (nframes, nod2D)
        surf2 = data2[name][:, :, 0]
        diff = surf1 - surf2

        # Shared color range for the two field panels
        vmin = min(surf1.min(), surf2.min())
        vmax = max(surf1.max(), surf2.max())
        if vmin == vmax:
            vmin -= 0.5
            vmax += 0.5

        # Symmetric color range for difference
        dmax = max(abs(diff.min()), abs(diff.max()))
        if dmax == 0:
            dmax = 1e-10

        frames_for_gif = []
        for step in range(nframes):
            fig, axes = plt.subplots(1, 3, figsize=(22, 6))

            # Panel 1: dir1
            tpc1 = axes[0].tripcolor(triang, surf1[step], shading="flat",
                                     cmap="RdYlBu_r", vmin=vmin, vmax=vmax)
            fig.colorbar(tpc1, ax=axes[0], shrink=0.8)
            axes[0].set_title(f"{name.capitalize()} [{label1}]")
            axes[0].set_aspect("equal")
            axes[0].set_xlabel("x (m)")
            axes[0].set_ylabel("y (m)")

            # Panel 2: dir2
            tpc2 = axes[1].tripcolor(triang, surf2[step], shading="flat",
                                     cmap="RdYlBu_r", vmin=vmin, vmax=vmax)
            fig.colorbar(tpc2, ax=axes[1], shrink=0.8)
            axes[1].set_title(f"{name.capitalize()} [{label2}]")
            axes[1].set_aspect("equal")
            axes[1].set_xlabel("x (m)")

            # Panel 3: difference
            tpc3 = axes[2].tripcolor(triang, diff[step], shading="flat",
                                     cmap="RdBu_r", vmin=-dmax, vmax=dmax)
            fig.colorbar(tpc3, ax=axes[2], shrink=0.8)
            axes[2].set_title(f"Diff ({label1} - {label2})")
            axes[2].set_aspect("equal")
            axes[2].set_xlabel("x (m)")

            fig.suptitle(f"Step {step}/{nsteps}", fontsize=14, fontweight="bold")
            fig.tight_layout()

            fname = os.path.join(img_dir, f"{name}_compare_step{step:04d}.png")
            fig.savefig(fname, dpi=100, bbox_inches="tight")
            frames_for_gif.append(fname)
            plt.close(fig)

        print(f"  Saved {nframes} comparison frames for '{name}' to {img_dir}/")

        # Print per-step difference statistics
        print(f"  {name} difference stats ({label1} - {label2}):")
        for step in range(nframes):
            d = diff[step]
            print(f"    step {step:3d}: max|diff|={np.abs(d).max():.6e}, "
                  f"mean|diff|={np.abs(d).mean():.6e}, rms={np.sqrt((d**2).mean()):.6e}")

        try:
            from PIL import Image
            imgs = [Image.open(f) for f in frames_for_gif]
            gif_path = os.path.join(output_dir, f"{name}_compare.gif")
            imgs[0].save(gif_path, save_all=True, append_images=imgs[1:],
                         duration=300, loop=0)
            print(f"  Animation saved to {gif_path}")
        except ImportError:
            print("  (Pillow not installed — skipping GIF animation)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot FESOM2 tracer scalar fields or compare two runs"
    )
    parser.add_argument("dirs", nargs="*", default=["mesh_output"],
                        help="Output directory (or two directories with --compare)")
    parser.add_argument("--compare", action="store_true",
                        help="Compare two directories side by side with difference")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory for comparison plots")
    args = parser.parse_args()

    if args.compare:
        if len(args.dirs) != 2:
            print("Error: --compare requires exactly two directories", file=sys.stderr)
            sys.exit(1)
        for d in args.dirs:
            if not os.path.isdir(d):
                print(f"Error: directory '{d}' not found", file=sys.stderr)
                sys.exit(1)
            for required in ["mesh_info.txt", "scalar_info.txt"]:
                if not os.path.isfile(os.path.join(d, required)):
                    print(f"Error: '{required}' not found in '{d}'", file=sys.stderr)
                    sys.exit(1)

        for d in args.dirs:
            info, data = load_scalars(d)
            wp_label = get_precision_label(info)
            print(f"Loaded scalars from '{d}' [{wp_label}]:")
            print(f"  Nodes: {info['nod2D']}, Layers: {info['nz']}, "
                  f"Steps: {info['nsteps']}, Tracers: {', '.join(info['tracer_names'])}")
            for name in info["tracer_names"]:
                arr = data[name]
                print(f"    {name}: min={arr.min():.6f}, max={arr.max():.6f}")

        print("\nComparing...")
        plot_compare(args.dirs[0], args.dirs[1], output_dir=args.output_dir)
        print("Done.")
    else:
        mesh_dir = args.dirs[0]
        if not os.path.isdir(mesh_dir):
            print(f"Error: directory '{mesh_dir}' not found", file=sys.stderr)
            sys.exit(1)
        for required in ["mesh_info.txt", "scalar_info.txt"]:
            if not os.path.isfile(os.path.join(mesh_dir, required)):
                print(f"Error: '{required}' not found in '{mesh_dir}'", file=sys.stderr)
                print("Make sure to run with both --save-mesh and --save-scalars", file=sys.stderr)
                sys.exit(1)

        info, data = load_scalars(mesh_dir)
        wp_label = get_precision_label(info)
        print(f"Loaded scalars from '{mesh_dir}' [{wp_label}]:")
        print(f"  Nodes:    {info['nod2D']}")
        print(f"  Layers:   {info['nz']}")
        print(f"  Steps:    {info['nsteps']} (+1 initial = {info['nsteps']+1} frames)")
        print(f"  Tracers:  {info['num_tracers']} ({', '.join(info['tracer_names'])})")
        for name in info["tracer_names"]:
            arr = data[name]
            print(f"    {name}: min={arr.min():.6f}, max={arr.max():.6f}")

        print("\nPlotting...")
        plot_scalars(mesh_dir)
        print("Done.")
