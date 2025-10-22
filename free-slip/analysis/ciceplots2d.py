#!/usr/bin/env python3

import sys
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature

if len(sys.argv) != 6:
    print("ciceplots2d_cartopy.py requires 5 arguments")
    print('  1. field name in file, ie. "aice"')
    print('  2. cice history file full path, ie. "/path/to/iceh.2012-03.nc"')
    print('  3. case name, used to annotate plot, ie. "CICE6.5.1"')
    print('  4. notes, used to annotate plot, ie. 2012 "March Mean"')
    print('  5. file string, used to create unique png filenames, ie. "Mar12"')
    sys.exit(1)

field = sys.argv[1]
pathf = sys.argv[2]
casen = sys.argv[3]
notes = sys.argv[4]
fstr  = sys.argv[5]
fname = os.path.basename(pathf)
title = field + " " + notes
cfnam = casen + " " + fname

# --- Read data ---
ds = Dataset(pathf, "r")

def getvar(name):
    if name in ds.variables:
        return ds.variables[name][:]
    # case-insensitive fallback
    for k in ds.variables:
        if k.lower() == name.lower():
            return ds.variables[k][:]
    raise KeyError(f"Variable '{name}' not found")

lons = getvar("TLON")
lats = getvar("TLAT")

var = getvar(field)
# expect [time,y,x], pick first time if 3D
if var.ndim == 3:
    var = var[0, :, :]
elif var.ndim == 2:
    var = var[:, :]
else:
    raise ValueError(f"Unsupported variable dims {var.shape} for '{field}'")

# convert zeros to NaN (to match original behavior)
var = np.array(var, dtype=float)
var[var == 0.0] = np.nan

# --- Colormap / normalization ---
cmap = plt.get_cmap("jet")
barticks = None
norm = None

if field in ["hi"]:
    bounds = np.arange(0, 2.05, 0.1)
    bounds = np.append(bounds, [2.25,2.5,2.75,3.0,3.25,3.5,3.75,4.0])
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend="max")
    barticks = [0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0]

if field in ["hs"]:
    bounds = np.arange(0, 1.02, 0.05)
    bounds = np.append(bounds, [1.5,2.0,2.5,3.0,3.5,4.0])
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend="max")
    barticks = [0,0.25,0.5,0.75,1.0,2.0,3.0,4.0]

def add_common(ax, title, cf_label):
    # Land/ocean styling similar to Basemap fillcontinents black / map boundary white
    ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=0)
    ax.add_feature(cfeature.LAND, facecolor="black", zorder=1)
    ax.coastlines(resolution="110m", linewidth=0.3, zorder=2)

    # Gridlines (labels as in original: only on global plot; for polars we'll skip labels)
    gl = ax.gridlines(draw_labels=True, linewidth=0.4, linestyle="--", color="0.5")
    try:
        gl.top_labels = False
        gl.right_labels = False
    except Exception:
        pass

    # Colorbar label set by caller (after plotting)
    ax.set_title(title)

def draw_scalar(ax, projection_name):
    # We are plotting lon/lat data; tell Cartopy the source CRS
    transform = ccrs.PlateCarree()
    # Use pcolormesh-like scatter with tiny points to mimic original; but pcolormesh is faster if gridded.
    # Try pcolormesh first; fall back to scatter if masked shapes cause issues.
    try:
        cf = ax.pcolormesh(lons, lats, var, transform=transform, cmap=cmap, shading="auto", norm=norm)
    except Exception:
        cf = ax.scatter(lons, lats, c=var, s=0.2, cmap=cmap, norm=norm, transform=transform)

    cbar = plt.colorbar(cf, ax=ax, shrink=0.8, pad=0.04)
    if barticks is not None:
        cbar.set_ticks(barticks)
    cbar.set_label(field)

    # annotate (same placement as original)
    ax.text(x=0.0, y=-0.1 if projection_name=="global" else -0.02,
            s=cfnam, transform=ax.transAxes,
            ha="left", va="top", fontsize="x-small")

# --- Figure 1: Global (cylindrical) ---
fig = plt.figure(figsize=(6,4), dpi=300)
ax = plt.axes(projection=ccrs.PlateCarree())
add_common(ax, title, field)
# extent 0..360 is fine; Cartopy expects degrees.
ax.set_global()
draw_scalar(ax, projection_name="global")
oname = f"figs/{field}_gl_{fstr}.png"
print("Saving file to", oname)
plt.savefig(oname, bbox_inches="tight")
plt.close(fig)

# --- Figure 2: North Polar Stereographic ---
fig = plt.figure(figsize=(6,4), dpi=300)
ax = plt.axes(projection=ccrs.NorthPolarStereo(central_longitude=-45))
# Limit to 45N..90N
ax.set_extent([-180, 180, 45, 90], crs=ccrs.PlateCarree())
# Turn off gridline labels on polar plots
gl = ax.gridlines(draw_labels=False, linewidth=0.4, linestyle="--", color="0.5")
ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=0)
ax.add_feature(cfeature.LAND, facecolor="black", zorder=1)
ax.coastlines(resolution="110m", linewidth=0.3, zorder=2)
ax.set_title(title)
draw_scalar(ax, projection_name="nh")
oname = f"figs/{field}_nh_{fstr}.png"
print("Saving file to", oname)
plt.savefig(oname, bbox_inches="tight")
plt.close(fig)

# --- Figure 3: South Polar Stereographic ---
fig = plt.figure(figsize=(6,4), dpi=300)
ax = plt.axes(projection=ccrs.SouthPolarStereo(central_longitude=180))
# Limit to 45S..90S
ax.set_extent([-180, 180, -90, -45], crs=ccrs.PlateCarree())
gl = ax.gridlines(draw_labels=False, linewidth=0.4, linestyle="--", color="0.5")
ax.add_feature(cfeature.OCEAN, facecolor="white", zorder=0)
ax.add_feature(cfeature.LAND, facecolor="black", zorder=1)
ax.coastlines(resolution="110m", linewidth=0.3, zorder=2)
ax.set_title(title)
draw_scalar(ax, projection_name="sh")
oname = f"figs/{field}_sh_{fstr}.png"
print("Saving file to", oname)
plt.savefig(oname, bbox_inches="tight")
plt.close(fig)

ds.close()
