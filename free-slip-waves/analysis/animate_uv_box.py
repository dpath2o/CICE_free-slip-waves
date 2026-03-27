# This script generates two MP4 animations (u and v on the U-grid) with a fixed colorbar
# range [-0.1, 0.1]. It is NaN-aware and reads a *range* of instantaneous files
# (iceh_inst.YYYY-MM-DD-hhmmss.nc) between --start and --end (inclusive).
#
# Usage example (run this on Gadi where your files live):
#   python make_uv_movies_fixed_scale.py \
#     --hist /g/data/gv90/da1339/cice-dirs/runs/RESULTS/free-slip/history \
#     --start 2005-01-01 --end 2005-01-03 --label free
#
# It will write: free_u_fixed.mp4 and free_v_fixed.mp4 in the current folder.
#
# Notes:
# - Only uses matplotlib (no seaborn).
# - One figure per movie (no subplots).
# - Colorbar is fixed to [-0.1, 0.1] as requested.
# - If ffmpeg is not on PATH, the script will fall back to Pillow writer (slower, larger).
#
#!/usr/bin/env python3
import os
import re
import glob
import argparse
from datetime import datetime

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from netCDF4 import Dataset

FNAME_RE = re.compile(r"iceh_inst\.(\d{4}-\d{2}-\d{2})-(\d{5,6})\.nc$")

def collect_files_in_range(hist_dir, start_date, end_date):
    pats = [
        os.path.join(hist_dir, "iceh_inst.*.nc"),
        os.path.join(hist_dir, "*", "iceh_inst.*.nc"),
    ]
    files = []
    for pat in pats:
        files.extend(glob.glob(pat))
    selected = []
    for f in files:
        m = FNAME_RE.search(os.path.basename(f))
        if not m:
            continue
        dstr, tstr = m.groups()
        if len(tstr) == 5:
            tstr = tstr.zfill(6)
        try:
            ts = datetime.strptime(f"{dstr}-{tstr}", "%Y-%m-%d-%H%M%S")
        except Exception:
            continue
        if start_date <= ts.date() <= end_date:
            selected.append((ts, f))
    selected.sort(key=lambda x: x[0])
    return [f for _, f in selected]

def read_uv_u(nc_path):
    with Dataset(nc_path, "r") as ds:
        def get(name_list):
            for n in name_list:
                if n in ds.variables:
                    return ds.variables[n][:]
            return None
        def sq(a):
            a = np.asarray(a)
            if a.ndim == 3:
                a = a[0]
            return a.astype(float)
        u = get(["uvel_1", "uvelU_1", "UVEL_U", "uvel"])
        v = get(["vvel_1", "vvelU_1", "VVEL_U", "vvel"])
        if u is None or v is None:
            raise RuntimeError(f"Missing uvel_1/vvel_1 in {os.path.basename(nc_path)}")
        return sq(u), sq(v)

def make_movie(frames_array, out_path, label, comp_label, vmin=-0.1, vmax=0.1, fps=10):
    fig = plt.figure(figsize=(6, 5), dpi=140)
    ax = fig.add_axes([0.10, 0.10, 0.78, 0.80])
    first = np.ma.masked_invalid(frames_array[0])
    im = ax.imshow(first, origin="lower", vmin=vmin, vmax=vmax, interpolation="nearest")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f"{comp_label} (m/s)")
    ax.set_title(f"{label}: {comp_label} on U-grid")
    ax.set_xlabel("i")
    ax.set_ylabel("j")

    def update(k):
        im.set_data(np.ma.masked_invalid(frames_array[k]))
        return (im,)

    try:
        writer = animation.FFMpegWriter(fps=fps, bitrate=2400)
        use_writer = writer
    except Exception:
        use_writer = animation.PillowWriter(fps=fps)
    ani = animation.FuncAnimation(fig, update, frames=len(frames_array), blit=True)
    ani.save(out_path, writer=use_writer)
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description="Make fixed-scale animations of uvel_1 and vvel_1 on the U-grid.")
    ap.add_argument("--hist", required=True, help="history directory containing iceh_inst.*.nc files")
    ap.add_argument("--start", required=True, help="start date YYYY-MM-DD (inclusive)")
    ap.add_argument("--end", required=True, help="end date YYYY-MM-DD (inclusive)")
    ap.add_argument("--label", default="exp", help="short label used in titles and output filenames")
    ap.add_argument("--vmin", type=float, default=-0.1, help="colorbar min (default -0.1)")
    ap.add_argument("--vmax", type=float, default=0.1, help="colorbar max (default 0.1)")
    ap.add_argument("--fps", type=int, default=10, help="frames per second (default 10)")
    args = ap.parse_args()

    start_date = datetime.strptime(args.start, "%Y-%m-%d").date()
    end_date   = datetime.strptime(args.end,   "%Y-%m-%d").date()

    files = collect_files_in_range(args.hist, start_date, end_date)
    if not files:
        print(f"No iceh_inst files found in range {args.start}..{args.end} under {args.hist}")
        return

    u_frames, v_frames = [], []
    for f in files:
        try:
            u, v = read_uv_u(f)
        except Exception as e:
            print(f"Skipping {os.path.basename(f)}: {e}")
            continue
        u_frames.append(u)
        v_frames.append(v)

    if not u_frames:
        print("No usable files/variables to animate.")
        return

    out_u = f"{args.label}_u_fixed.mp4"
    out_v = f"{args.label}_v_fixed.mp4"
    make_movie(u_frames, out_u, args.label, "u", vmin=args.vmin, vmax=args.vmax, fps=args.fps)
    make_movie(v_frames, out_v, args.label, "v", vmin=args.vmin, vmax=args.vmax, fps=args.fps)
    print(f"Wrote {out_u} and {out_v}")

if __name__ == "__main__":
    main()
