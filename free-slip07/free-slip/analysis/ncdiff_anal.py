#!/usr/bin/env python3
"""
ncdiff_anal.py

Compute ncdiff (A - B by default) and analyse the resulting NetCDF difference file.

Example:
  python3 ncdiff_anal.py /path/to/A.nc /path/to/B.nc /path/to/diff.nc --summary-dir /some/dir --top 6

Notes:
- Requires NCO's `ncdiff` in your PATH.
- Analysis is NaN-aware; stats are computed on finite values only.
- Edge stats assume dims named nj (y) and ni (x). If absent, edge stats are skipped.

Author: dpath2o, Oct 2025
"""
import argparse
import os
import sys
import subprocess as sp
from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

# -----------------------
# Helpers
# -----------------------
def run_ncdiff(A, B, diff_out, reverse=False, overwrite=False):
    if not shutil_which("ncdiff"):
        sys.exit("ERROR: `ncdiff` not found in PATH. Load NCO module or install NCO.")
    if Path(diff_out).exists() and not overwrite:
        print(f"[ncdiff] {diff_out} already exists (use --overwrite to re-create). Skipping ncdiff.")
        return
    cmd = ["ncdiff", "-O", B, A, diff_out] if reverse else ["ncdiff", "-O", A, B, diff_out]
    print("[ncdiff] Running:", " ".join(cmd))
    res = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
    if res.returncode != 0:
        print(res.stdout)
        print(res.stderr, file=sys.stderr)
        sys.exit(f"ERROR: ncdiff failed with code {res.returncode}")
    print(f"[ncdiff] Wrote: {diff_out}")

def shutil_which(name: str):
    from shutil import which
    return which(name)

def open_ds(path):
    """Open a dataset, falling back to netcdf4 engine if necessary."""
    try:
        return xr.open_dataset(path, decode_times=False)
    except Exception:
        return xr.open_dataset(path, engine="netcdf4", decode_times=False)

def diff_stats(var_da, time_dim="time", j_dim="nj", i_dim="ni"):
    arr         = var_da.values
    total_count = int(arr.size) if arr.size else 0
    if total_count == 0:
        return dict(var=var_da.name, count=0, nonzero=0,
                    max_abs=np.nan, mean_abs=np.nan, rmse=np.nan,
                    max_time_index=np.nan, max_j=np.nan, max_i=np.nan,
                    edge_frac=np.nan, edge_count=0)
    finite = np.isfinite(arr)
    if not finite.any():
        return dict(var=var_da.name, count=total_count, nonzero=0,
                    max_abs=np.nan, mean_abs=np.nan, rmse=np.nan,
                    max_time_index=np.nan, max_j=np.nan, max_i=np.nan,
                    edge_frac=np.nan, edge_count=0)
    abs_arr  = np.abs(arr)
    nz_mask  = finite & (abs_arr > 0)
    nz_count = int(nz_mask.sum())
    if nz_count == 0:
        return dict(var=var_da.name, count=total_count, nonzero=0,
                    max_abs=0.0, mean_abs=0.0, rmse=0.0,
                    max_time_index=np.nan, max_j=np.nan, max_i=np.nan,
                    edge_frac=0.0, edge_count=0)

    nz_vals  = abs_arr[nz_mask]
    max_abs  = float(np.nanmax(nz_vals))
    mean_abs = float(np.nanmean(nz_vals))
    rmse     = float(np.sqrt(np.nanmean(nz_vals**2)))
    # max location
    max_pos        = np.unravel_index(np.nanargmax(abs_arr), arr.shape)
    coords         = list(var_da.dims)
    max_time_index = np.nan; max_j = np.nan; max_i = np.nan
    for ax, dim in enumerate(coords):
        idx = max_pos[ax]
        if dim == time_dim:
            max_time_index = int(idx)
        elif dim == j_dim:
            max_j = int(idx)
        elif dim == i_dim:
            max_i = int(idx)
    # edge fraction
    edge_frac = np.nan; edge_count = 0
    if j_dim in coords and i_dim in coords:
        j_axis    = coords.index(j_dim)
        i_axis    = coords.index(i_dim)
        shape     = arr.shape
        edge_mask = np.zeros_like(arr, dtype=bool)
        slicer    = [slice(None)] * arr.ndim
        s = slicer.copy(); s[j_axis] = 0;              edge_mask[tuple(s)] = True
        s = slicer.copy(); s[j_axis] = shape[j_axis]-1; edge_mask[tuple(s)] = True
        s = slicer.copy(); s[i_axis] = 0;              edge_mask[tuple(s)] = True
        s = slicer.copy(); s[i_axis] = shape[i_axis]-1; edge_mask[tuple(s)] = True
        edge_count = int((edge_mask & nz_mask).sum())
        edge_frac  = float(edge_count / nz_count) if nz_count else 0.0
    return dict(var=var_da.name, count=total_count, nonzero=nz_count,
                max_abs=max_abs, mean_abs=mean_abs, rmse=rmse,
                max_time_index=max_time_index, max_j=max_j, max_i=max_i,
                edge_frac=edge_frac, edge_count=edge_count)

def edge_maxima(var_da, j_dim="nj", i_dim="ni"):
    """
    Return per-edge maxima (abs) if both j_dim and i_dim exist.
    Keys: west_i0, east_i-1, south_j0, north_j-1
    """
    out = {}
    if j_dim not in var_da.dims or i_dim not in var_da.dims:
        return out
    arr    = np.abs(var_da.values)
    finite = np.isfinite(arr)
    if not finite.any():
        return out
    # Build slicers
    dims  = list(var_da.dims)
    j_ax  = dims.index(j_dim)
    i_ax  = dims.index(i_dim)
    shape = arr.shape
    def max_on_slice(jslice=None, islice=None):
        sl = [slice(None)] * arr.ndim
        if jslice is not None: sl[j_ax] = jslice
        if islice is not None: sl[i_ax] = islice
        sub = arr[tuple(sl)]
        sub = sub[np.isfinite(sub)]
        return float(np.nanmax(sub)) if sub.size else np.nan
    out["west_i0"]   = max_on_slice(islice=0)
    out["east_i-1"]  = max_on_slice(islice=shape[i_ax]-1)
    out["south_j0"]  = max_on_slice(jslice=0)
    out["north_j-1"] = max_on_slice(jslice=shape[j_ax]-1)
    return out

def select_for_plot(da, row, time_dim="time", j_dim="nj", i_dim="ni"):
    """Slice to the time index of max (if any) and first index of extra dims, for plotting."""
    idx = {}
    if time_dim in da.dims and np.isfinite(row.get("max_time_index", np.nan)):
        idx[time_dim] = int(row["max_time_index"])
    for d in da.dims:
        if d not in (time_dim, j_dim, i_dim):
            idx[d] = 0
    return da.isel(**idx) if idx else da

def make_plot(da2, title, out_png):
    arr = np.abs(da2.values)
    vmax = float(np.nanmax(arr)) if np.isfinite(arr).any() else 1.0
    if not np.isfinite(vmax) or vmax == 0.0:
        vmax = 1.0
    fig = plt.figure(figsize=(6, 5))
    ax  = plt.gca()
    im  = ax.imshow(arr, origin="lower", vmin=0.0, vmax=vmax)
    ax.set_title(title)
    # Prefer dims if present
    j_dim = next((d for d in da2.dims if d in ("nj","y","j")), "y")
    i_dim = next((d for d in da2.dims if d in ("ni","x","i")), "x")
    ax.set_xlabel(i_dim); ax.set_ylabel(j_dim)
    cb = plt.colorbar(im); cb.set_label("|difference|")
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)

# -----------------------
# Main
# -----------------------
def main():
    p = argparse.ArgumentParser(description="Run ncdiff on two NetCDF files and analyse the diff (NaN-aware).")
    p.add_argument("A", help="Path to file A (minuend)")
    p.add_argument("B", help="Path to file B (subtrahend)")
    p.add_argument("DIFF", help="Output path for ncdiff result")
    p.add_argument("--reverse", action="store_true", help="Compute B - A instead of A - B")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing DIFF file")
    p.add_argument("--summary-dir", default=None, help="Directory for CSV/TXT summaries and plots (default: DIFF's directory)")
    p.add_argument("--top", type=int, default=6, help="Max number of variables to plot (with nonzero diffs)")
    p.add_argument("--edge-vars", nargs="*", help="Optional variable names to print extra per-edge max(|delta|) (requires dims nj,ni)")
    args = p.parse_args()

    A = Path(args.A).expanduser().resolve()
    B = Path(args.B).expanduser().resolve()
    DIFF = Path(args.DIFF).expanduser().resolve()
    if not A.is_file(): sys.exit(f"ERROR: A not found: {A}")
    if not B.is_file(): sys.exit(f"ERROR: B not found: {B}")
    out_dir = Path(args.summary_dir) if args.summary_dir else DIFF.parent
    out_dir = out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    # 1) ncdiff
    run_ncdiff(str(A), str(B), str(DIFF), reverse=args.reverse, overwrite=args.overwrite)
    # 2) Analysis
    print(f"[analysis] Opening diff file: {DIFF}")
    try:
        ds = open_ds(str(DIFF))
    except Exception as e:
        sys.exit(f"ERROR: failed to open diff file: {e}")
    dims = dict(ds.dims)
    time_dim = "time" if "time" in ds.dims else None
    j_dim = "nj" if "nj" in ds.dims else None
    i_dim = "ni" if "ni" in ds.dims else None
    # Collect numeric variables
    vars_to_check = [v for v in ds.data_vars if ds[v].dtype.kind in ("f","i","u")]
    if not vars_to_check:
        print("[analyze] No numeric data variables found. Exiting.")
        sys.exit(0)
    print(f"[analyze] Numeric variables to check ({len(vars_to_check)}): {', '.join(vars_to_check)}")
    # Stats per var
    records = []
    for v in vars_to_check:
        rec = diff_stats(ds[v], time_dim=time_dim or "time", j_dim=j_dim or "nj", i_dim=i_dim or "ni")
        records.append(rec)
    df = pd.DataFrame(records).sort_values(["nonzero","max_abs"], ascending=[False, False])
    # Write CSV + TXT summaries
    csv_path = out_dir / (DIFF.stem + "_diff_summary.csv")
    txt_path = out_dir / (DIFF.stem + "_diff_summary.txt")
    df.to_csv(csv_path, index=False)
    with open(txt_path, "w") as f:
        f.write("NetCDF diff analysis\n")
        f.write(f"Diff file: {DIFF}\n")
        f.write(f"Dims: {dims}\n")
        f.write(f"Checked variables: {vars_to_check}\n\n")
        f.write(df.to_string(index=False))
        f.write("\n")
    print(f"[analysis] Wrote: {csv_path}")
    print(f"[analysis] Wrote: {txt_path}")
    # Plots for top-N with differences
    plot_dir = out_dir / (DIFF.stem + "_plots")
    plot_dir.mkdir(parents=True, exist_ok=True)
    created = []
    top_vars = df[df["nonzero"] > 0].head(max(0, args.top))["var"].tolist()
    for v in top_vars:
        row = df[df["var"] == v].iloc[0].to_dict()
        da = ds[v]
        da2 = select_for_plot(da, row, time_dim=time_dim or "time", j_dim=j_dim or "nj", i_dim=i_dim or "ni")
        # If no 2D (j,i) slice possible, skip plot
        if not ((j_dim in da2.dims if j_dim else False) and (i_dim in da2.dims if i_dim else False)):
            continue
        t_str = f"{int(row['max_time_index'])}" if np.isfinite(row.get("max_time_index", np.nan)) else "all"
        out_png = plot_dir / f"diff_{v}_t{t_str}.png"
        title = f"|Δ| for {v} (t={t_str})"
        make_plot(da2, title, out_png)
        created.append(str(out_png))
    if created:
        print(f"[analysis] Plots: {len(created)} files in {plot_dir}")
    else:
        print("[analysis] No plots created (no variables with non-zero diffs or no 2D (nj,ni) slice).")
    # Optional edge sanity check
    if args.edge_vars:
        edge_txt = out_dir / (DIFF.stem + "_edge_max.txt")
        with open(edge_txt, "w") as f:
            f.write(f"Edge maxima (abs diffs) for selected variables\n")
            f.write(f"File: {DIFF}\n\n")
            for v in args.edge_vars:
                if v not in ds.variables:
                    f.write(f"{v}: (not found)\n")
                    continue
                em = edge_maxima(ds[v], j_dim=j_dim or "nj", i_dim=i_dim or "ni")
                if not em:
                    f.write(f"{v}: (edge stats not available; dims missing)\n")
                else:
                    f.write(f"{v}: west_i0={em['west_i0']:.6e}, "
                            f"east_i-1={em['east_i-1']:.6e}, "
                            f"south_j0={em['south_j0']:.6e}, "
                            f"north_j-1={em['north_j-1']:.6e}\n")
        print(f"[analysis] Wrote edge check: {edge_txt}")
    print("[done]")

if __name__ == "__main__":
    main()
