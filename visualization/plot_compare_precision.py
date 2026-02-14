#!/usr/bin/env python3
"""Compare FESOM2 tracer scalar fields across multiple precision variants.

Usage:
    python plot_compare_precision.py                    # auto-detect build_nvidia_*/mesh_output
    python plot_compare_precision.py dir1 dir2 dir3     # explicit directories
    python plot_compare_precision.py --diff dir1 dir2 dir3  # add difference panels (first dir = reference)
    python plot_compare_precision.py --diff --ref-label DP build_nvidia_dp/mesh_output ...

Plots surface (k=0) tracer fields side by side for each timestep.
With --diff, adds difference panels using the first directory as reference.
"""

import sys
import os
import argparse
import numpy as np


def load_mesh(dirpath):
    info = {}
    with open(os.path.join(dirpath, "mesh_info.txt")) as f:
        for line in f:
            key, val = line.strip().split("=")
            info[key] = int(val)

    nod2D = info["nod2D"]
    elem2D = info["elem2D"]

    coords = np.fromfile(
        os.path.join(dirpath, "coord_nod2D.bin"), dtype=np.float64
    ).reshape((2, nod2D), order="F")
    x = coords[0, :]
    y = coords[1, :]

    triangles = np.fromfile(
        os.path.join(dirpath, "elem2D_nodes.bin"), dtype=np.int32
    ).reshape((3, elem2D), order="F").T

    return x, y, triangles


def load_scalars(dirpath):
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
        data[name] = raw.reshape((nframes, nod2D, nz))

    return info, data


def get_precision_label(info):
    wp = info.get("wp", 8)
    return {2: "HP", 4: "SP", 8: "DP"}.get(wp, f"WP{wp}")


def make_triangulation(x, y, triangles):
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

    return triang


def plot_multi_precision(dirs, labels, show_diff, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ndirs = len(dirs)

    # Load all data
    meshes = [load_mesh(d) for d in dirs]
    scalar_sets = [load_scalars(d) for d in dirs]
    infos = [s[0] for s in scalar_sets]
    datas = [s[1] for s in scalar_sets]

    # Auto-labels from precision if not provided
    if labels is None:
        labels = [get_precision_label(info) for info in infos]

    # Validate compatible grids
    ref_nod = infos[0]["nod2D"]
    ref_steps = infos[0]["nsteps"]
    for i, info in enumerate(infos[1:], 1):
        if info["nod2D"] != ref_nod:
            raise ValueError(f"Node count mismatch: dir0={ref_nod} vs dir{i}={info['nod2D']}")
        if info["nsteps"] != ref_steps:
            raise ValueError(f"Step count mismatch: dir0={ref_steps} vs dir{i}={info['nsteps']}")

    nsteps = ref_steps
    nframes = nsteps + 1
    tracer_names = infos[0]["tracer_names"]

    # Use first mesh for triangulation (same grid)
    triang = make_triangulation(*meshes[0])

    os.makedirs(output_dir, exist_ok=True)
    img_dir = os.path.join(output_dir, "plots")
    os.makedirs(img_dir, exist_ok=True)

    # Layout: with --diff, 2 rows (fields + diffs); without, 1 row
    ncols = ndirs
    nrows = 2 if show_diff else 1
    for name in tracer_names:
        surfaces = [datas[i][name][:, :, 0] for i in range(ndirs)]  # (nframes, nod2D)

        # Shared color range across all precision variants
        vmin = min(s.min() for s in surfaces)
        vmax = max(s.max() for s in surfaces)
        if vmin == vmax:
            vmin -= 0.5
            vmax += 0.5

        # Difference arrays (relative to first = reference)
        diffs = []
        dmax = 1e-15
        if show_diff:
            for i in range(1, ndirs):
                diffs.append(surfaces[i] - surfaces[0])
            dmax = max(max(np.abs(d).max() for d in diffs), dmax)

        frames_for_gif = []
        for step in range(nframes):
            fig_w = 6 * ncols + 1
            fig_h = 5.5 * nrows + 0.8
            fig, axes = plt.subplots(nrows, ncols, figsize=(fig_w, fig_h),
                                     squeeze=False)

            # Row 0: field panels
            for i in range(ndirs):
                ax = axes[0][i]
                tpc = ax.tripcolor(triang, surfaces[i][step], shading="flat",
                                   cmap="RdYlBu_r", vmin=vmin, vmax=vmax)
                fig.colorbar(tpc, ax=ax, shrink=0.8)
                ax.set_title(f"{name.capitalize()} [{labels[i]}]", fontsize=12)
                ax.set_aspect("equal")
                ax.set_xlabel("x (m)")
                if i == 0:
                    ax.set_ylabel("y (m)")

            # Row 1: difference panels (skip col 0 = reference)
            if show_diff:
                # Leave first column of row 1 empty (reference has no diff)
                axes[1][0].set_visible(False)
                for j, d in enumerate(diffs):
                    ax = axes[1][j + 1]
                    tpc = ax.tripcolor(triang, d[step], shading="flat",
                                       cmap="RdBu_r", vmin=-dmax, vmax=dmax)
                    fig.colorbar(tpc, ax=ax, shrink=0.8, format="%.2e")
                    ax.set_title(f"{labels[j+1]} - {labels[0]}", fontsize=12)
                    ax.set_aspect("equal")
                    ax.set_xlabel("x (m)")

            fig.suptitle(f"Step {step}/{nsteps}", fontsize=14, fontweight="bold")
            fig.tight_layout()

            fname = os.path.join(img_dir, f"{name}_step{step:04d}.png")
            fig.savefig(fname, dpi=100, bbox_inches="tight")
            frames_for_gif.append(fname)
            plt.close(fig)

        print(f"  Saved {nframes} frames for '{name}' to {img_dir}/")

        # Print difference statistics
        if show_diff:
            for j, d in enumerate(diffs):
                print(f"\n  {name} diff stats ({labels[j+1]} - {labels[0]}):")
                for step in range(nframes):
                    ds = d[step]
                    print(f"    step {step:3d}: max|diff|={np.abs(ds).max():.6e}, "
                          f"rms={np.sqrt((ds**2).mean()):.6e}")

        # GIF animation
        try:
            from PIL import Image
            imgs = [Image.open(f) for f in frames_for_gif]
            gif_path = os.path.join(output_dir, f"{name}_compare.gif")
            imgs[0].save(gif_path, save_all=True, append_images=imgs[1:],
                         duration=400, loop=0)
            print(f"  Animation: {gif_path}")
        except ImportError:
            print("  (Pillow not installed — skipping GIF)")


def find_build_dirs(base):
    """Auto-detect build_nvidia_*/mesh_output directories, sorted DP > SP > HP."""
    candidates = []
    for name in sorted(os.listdir(base)):
        mesh_out = os.path.join(base, name, "mesh_output")
        if name.startswith("build_nvidia_") and os.path.isdir(mesh_out):
            scalar_info = os.path.join(mesh_out, "scalar_info.txt")
            if os.path.isfile(scalar_info):
                candidates.append(mesh_out)

    # Sort by wp descending (DP=8 first, then SP=4, then HP=2)
    def sort_key(d):
        with open(os.path.join(d, "scalar_info.txt")) as f:
            for line in f:
                if line.startswith("wp="):
                    return -int(line.strip().split("=")[1])
        return 0

    candidates.sort(key=sort_key)
    return candidates


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compare FESOM2 tracer scalars across precision variants"
    )
    parser.add_argument("dirs", nargs="*",
                        help="Output directories (mesh_output). If omitted, auto-detects build_nvidia_*/mesh_output")
    parser.add_argument("--diff", action="store_true",
                        help="Show difference panels (first dir = reference)")
    parser.add_argument("--labels", nargs="*", default=None,
                        help="Custom labels for each directory (default: auto from wp)")
    parser.add_argument("--output-dir", "-o", default="precision_comparison",
                        help="Output directory for plots (default: precision_comparison)")
    args = parser.parse_args()

    base = os.path.dirname(os.path.abspath(__file__))

    if not args.dirs:
        args.dirs = find_build_dirs(base)
        if not args.dirs:
            print("Error: no build_nvidia_*/mesh_output directories found", file=sys.stderr)
            sys.exit(1)
        print(f"Auto-detected {len(args.dirs)} output directories:")
        for d in args.dirs:
            print(f"  {d}")

    if len(args.dirs) < 2:
        print("Error: need at least 2 directories to compare", file=sys.stderr)
        sys.exit(1)

    # Validate
    for d in args.dirs:
        if not os.path.isdir(d):
            print(f"Error: '{d}' not found", file=sys.stderr)
            sys.exit(1)
        for req in ["mesh_info.txt", "scalar_info.txt"]:
            if not os.path.isfile(os.path.join(d, req)):
                print(f"Error: '{req}' not found in '{d}'", file=sys.stderr)
                sys.exit(1)

    # Load and summarize
    for d in args.dirs:
        info, data = load_scalars(d)
        label = get_precision_label(info)
        print(f"\n  [{label}] {d}")
        print(f"    Nodes: {info['nod2D']}, Layers: {info['nz']}, Steps: {info['nsteps']}")
        for name in info["tracer_names"]:
            arr = data[name]
            print(f"    {name}: min={arr.min():.6f}, max={arr.max():.6f}")

    print(f"\nPlotting to {args.output_dir}/...")
    plot_multi_precision(args.dirs, args.labels, args.diff, args.output_dir)
    print("Done.")
