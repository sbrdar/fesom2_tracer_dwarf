# Visualization

Python scripts for plotting mesh geometry, tracer fields, and multi-precision comparisons. All scripts live in `visualization/`.

## Environment setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install numpy matplotlib
pip install Pillow   # optional, for animated GIF export
```

## Generating output data

All plotting scripts read binary files written by the dwarf executables. Run with `--save-mesh` and/or `--save-scalars` to produce them:

```bash
cd build_gnu_dp
./run.sh 1 30 30 10 --periodic --save-mesh --save-scalars
```

This writes to `mesh_output/` (override with `--output-dir DIR`):

| File | Contents |
|---|---|
| `mesh_info.txt` | Node/element/level counts |
| `coord_nod2D.bin` | Node coordinates (float64, Fortran column-major) |
| `elem2D_nodes.bin` | Triangle connectivity (int32, 0-based) |
| `zbar.bin` | Vertical level depths (float64) |
| `depth.bin` | Bathymetry per node (float64) |
| `scalar_info.txt` | Scalar metadata (nod2D, nz, nsteps, wp, tracer names) |
| `<tracer>.bin` | Tracer field per frame (float64, shape: nframes x nod2D x nz) |

## Scripts

### `plot_mesh.py` — surface triangulation

Plots the 2D triangular mesh and prints mesh statistics (node/element count, coordinate ranges, depth range).

```bash
python visualization/plot_mesh.py [mesh_output_dir]
```

- Default directory: `mesh_output`
- Saves `mesh_plot.png` inside the mesh directory
- Shows the plot interactively

**Output**: a single PNG showing the surface triangulation with node count and element count in the title.

### `plot_scalars.py` — tracer field evolution

Plots the surface layer (k=0) of each tracer at every timestep. Supports both single-directory plotting and two-directory side-by-side comparison.

#### Single-directory mode

```bash
python visualization/plot_scalars.py [mesh_output_dir]
```

- Default directory: `mesh_output`
- Requires both `--save-mesh` and `--save-scalars` output in the directory
- Creates `scalar_plots/` subdirectory with per-step PNGs (`<tracer>_step0000.png`, ...)
- Each frame shows a tripcolor plot with a shared colorbar range across all steps
- Title includes precision label (DP/SP/HP) auto-detected from `scalar_info.txt`
- If Pillow is installed, saves `<tracer>_animation.gif` (300 ms per frame, looping)

#### Two-directory comparison mode

```bash
python visualization/plot_scalars.py --compare dir1 dir2 [--output-dir DIR]
```

- Requires exactly two directories, both with mesh and scalar data on the same grid
- Produces 3-panel figures per step: dir1 field, dir2 field, difference (dir1 - dir2)
- Field panels share a color range; difference panel uses a symmetric diverging colormap (`RdBu_r`)
- Prints per-step difference statistics: max|diff|, mean|diff|, RMS
- Saves frames to `compare_plots/` and optionally a `<tracer>_compare.gif`
- Default output directory: `compare_<label1>_vs_<label2>`

### `plot_compare_precision.py` — multi-precision comparison

Compares tracer fields across two or more precision variants side by side. Designed for DP vs SP vs HP comparisons.

```bash
# Auto-detect build_nvidia_*/mesh_output directories
python visualization/plot_compare_precision.py

# Explicit directories
python visualization/plot_compare_precision.py dir1 dir2 dir3

# With difference panels (first directory = reference)
python visualization/plot_compare_precision.py --diff dir1 dir2 dir3

# Custom labels and output directory
python visualization/plot_compare_precision.py --diff --labels DP SP HP -o my_comparison dir1 dir2 dir3
```

**Options**:

| Flag | Description |
|---|---|
| `--diff` | Add a second row of difference panels (each dir minus the first/reference) |
| `--labels L1 L2 ...` | Custom panel labels (default: auto-detected from `wp` in `scalar_info.txt`) |
| `-o`, `--output-dir` | Output directory (default: `precision_comparison`) |

**Auto-detection**: when invoked with no directory arguments, scans for `build_nvidia_*/mesh_output` directories and sorts them DP > SP > HP.

**Layout**:
- Without `--diff`: one row of N panels, one per directory
- With `--diff`: two rows — top row shows fields, bottom row shows difference from the reference (first directory). The reference column in the bottom row is left blank.

**Output**:
- `plots/` subdirectory with per-step PNGs
- Per-step difference statistics printed to stdout (max|diff|, RMS)
- Animated GIF per tracer if Pillow is available (400 ms per frame)

**Periodic mesh handling**: all three scripts automatically detect and mask wrap-around triangles on periodic meshes. Triangles spanning more than half the domain in x or y are masked to prevent visual artifacts.

## Full example workflow

```bash
# Build and run two precision variants
./configure.sh --compiler gnu --precision dp --clean --build
cd build_gnu_dp && ./run.sh 1 30 30 10 --periodic --save-mesh --save-scalars && cd ..

./configure.sh --compiler gnu --precision sp --clean --build
cd build_gnu_sp && ./run.sh 1 30 30 10 --periodic --save-mesh --save-scalars && cd ..

# Plot mesh
python visualization/plot_mesh.py build_gnu_dp/mesh_output

# Plot tracer evolution for one run
python visualization/plot_scalars.py build_gnu_dp/mesh_output

# Compare two runs side by side with difference
python visualization/plot_scalars.py --compare build_gnu_dp/mesh_output build_gnu_sp/mesh_output

# Multi-precision comparison with difference panels
python visualization/plot_compare_precision.py --diff \
    build_gnu_dp/mesh_output build_gnu_sp/mesh_output
```
