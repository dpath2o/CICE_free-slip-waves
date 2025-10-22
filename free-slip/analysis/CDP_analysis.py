# CDP quick-look: instantaneous u,v (faces) + daily mean KuxE,KuyN,uvel,vvel
# Works on the "BCtest" rectangular grid; robust to slightly different var names.

import os, glob, re
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib
matplotlib.use('Agg')   # must be set before importing pyplot
import matplotlib.pyplot as plt

# -------------------
# CONFIG
# -------------------
DATA_DIR   = os.path.expanduser('~/AFIM_runBox1/history')
outdir     = os.path.join(os.path.expanduser('~/AFIM_runBox1/'), 'CDP_analysis')
SIM_LABEL  = 'free-slip-no-CDP_ISOTROPIC_uniform-west'   # used in figure filenames
EDGE_FALLBACK_WIDTH = 2       # cells, used ONLY if F2/Kux/Kuy unavailable
EPS = 1e-20
rows = []

os.makedirs(outdir, exist_ok=True)

# -------------------
# HELPERS
# -------------------
def find_files(pattern):
    files = sorted(glob.glob(os.path.join(DATA_DIR, pattern)))
    if not files:
        raise FileNotFoundError(f'No files match {pattern} in {DATA_DIR}')
    return files

def pick_var(ds, preferred, variants=None):
    """Return the first present var among [preferred] + variants."""
    candidates = [preferred] + (variants or [])
    for k in ds.data_vars:
        for c in candidates:
            if k == c:
                return k
    # last chance: relaxed pattern match
    for k in ds.data_vars:
        if re.fullmatch(preferred+r'(_\d+)?', k):
            return k
    raise KeyError(f"None of {candidates} found in dataset.")

def dims_ij(da):
    """Return (nj,ni) dim names in order."""
    for cand_i in ('ni','xi','x','i'):
        for cand_j in ('nj','yi','y','j'):
            if cand_j in da.dims and cand_i in da.dims:
                return cand_j, cand_i
    # fallback: assume last two are (y,x)
    return da.dims[-2], da.dims[-1]

def load_inst():
    inst_files = find_files('iceh_inst.*.nc')
    # open lazily to keep memory modest
    ds = xr.open_mfdataset(inst_files, combine='by_coords', parallel=False)
    return ds

def load_daily():
    daily_files = find_files('iceh.2005-01-*.nc')
    ds = xr.open_mfdataset(daily_files, combine='by_coords', parallel=False)
    return ds

def build_edge_masks(daily):
    """Return boolean masks (E_mask, N_mask) marking where CDP is active on E and N faces."""
    # Prefer F2E/F2N if present
    try:
        f2e_name = pick_var(daily, 'F2E')
        F2E = daily[f2e_name].load()
    except Exception:
        F2E = None
    try:
        f2n_name = pick_var(daily, 'F2N')
        F2N = daily[f2n_name].load()
    except Exception:
        F2N = None

    # If not present, infer from |KuxE|/|KuyN|
    if F2E is None:
        try:
            kuxe_name = pick_var(daily, 'KuxE')
            F2E = xr.where(np.abs(daily[kuxe_name]) > EPS, 1.0, 0.0)
        except Exception:
            F2E = None
    if F2N is None:
        try:
            kuyn_name = pick_var(daily, 'KuyN')
            F2N = xr.where(np.abs(daily[kuyn_name]) > EPS, 1.0, 0.0)
        except Exception:
            F2N = None

    # Final fallback: synthetic top/bottom band for E, left/right band for N
    if F2E is None or F2N is None:
        # Create arrays using daily uvel as a template for shape
        templ = None
        for cand in ('uvel','vvel'):
            if cand in daily:
                templ = daily[cand].isel(time=0, drop=True)
                break
        if templ is None:
            # use any 2D var
            templ = next(iter(daily.data_vars.values())).isel(time=0, drop=True)
        nj_name, ni_name = dims_ij(templ)

        ny = templ.sizes[nj_name]
        nx = templ.sizes[ni_name]
        jj = xr.DataArray(np.arange(ny), dims=(nj_name,))
        ii = xr.DataArray(np.arange(nx), dims=(ni_name,))

        if F2E is None:
            band = ( (jj < EDGE_FALLBACK_WIDTH) | (jj >= ny-EDGE_FALLBACK_WIDTH) )
            F2E = xr.DataArray(band, dims=(nj_name,)).astype(float).broadcast_like(templ)
        if F2N is None:
            band = ( (ii < EDGE_FALLBACK_WIDTH) | (ii >= nx-EDGE_FALLBACK_WIDTH) )
            F2N = xr.DataArray(band, dims=(ni_name,)).astype(float).broadcast_like(templ)

    E_mask = (F2E.max('time') if 'time' in F2E.dims else F2E) > 0
    N_mask = (F2N.max('time') if 'time' in F2N.dims else F2N) > 0
    return E_mask, N_mask

def masked_mean_abs(da, mask):
    # broadcast mask to da
    m = mask.broadcast_like(da.isel(time=0, drop=True) if 'time' in da.dims else da)
    return float(np.abs(da).where(m, drop=True).mean().values)

def push(var, grid, edge_mean, int_mean, use_abs=False, note=''):
    rows.append(dict(
        sim=SIM_LABEL, var=var, grid=grid, abs=use_abs,
        edge_mean=edge_mean, interior_mean=int_mean,
        ratio=edge_mean/(int_mean+EPS), note=note
    ))

def bcast_mask(mask, target_da):
    """Broadcast boolean mask to match target_da's (time-less) shape."""
    base = target_da.isel(time=0, drop=True) if 'time' in target_da.dims else target_da
    return mask.broadcast_like(base)

def mean_over(mask, da, use_abs=False, time_mean=True):
    """Spatial mean over mask; optional |.| and time-mean first."""
    data = da
    if time_mean and 'time' in data.dims:
        data = data.mean('time')
    if use_abs:
        data = np.abs(data)
    M = bcast_mask(mask, data)
    return float(data.where(M).mean().values)

def hours_since_start(cf_times):
    t0 = cf_times[0]
    # cftime deltas behave like datetime.timedelta
    return np.array([ (t - t0).days*24.0 + (t - t0).seconds/3600.0 for t in cf_times ],dtype=float)

def numeric_time_seconds(da_time):
    th = hours_since_start(da_time)              # (hours)
    return th * 3600.0                           # seconds

def time_derivative(da):
    """Central finite-diff in time using real seconds, returns same dims."""
    tsec = numeric_time_seconds(da['time'].values)
    arr  = da.values
    dadt = np.gradient(arr, tsec, axis=0)
    return xr.DataArray(dadt, coords=da.coords, dims=da.dims,
                        name=f"d{da.name or 'var'}/dt")

def find_dx_dy_for_faces(template_da):
    """
    Best-effort grid spacing (meters) for E and N faces.
    If not found in files, fall back to 1.0 (index units).
    """
    cand_dx = [k for k in daily.data_vars if k.lower() in ('dx','dxu','dxe','dxt','hte','g_hte')]
    cand_dy = [k for k in daily.data_vars if k.lower() in ('dy','dyv','dyn','dyt','htn','g_htn')]
    DX = None; DY = None
    if cand_dx:
        try:
            DX = daily[cand_dx[0]].mean().item()
        except Exception:
            pass
    if cand_dy:
        try:
            DY = daily[cand_dy[0]].mean().item()
        except Exception:
            pass
    if DX is None: DX = 1.0
    if DY is None: DY = 1.0
    return float(DX), float(DY)

def get_spacing_line(ds, names, row_dim, col_dim, row_idx, prefer_row=True):
    """Return 1D spacing array along a line:
       - prefer_row=True  => take a row at fixed row_idx, vary along col_dim
       - prefer_row=False => take a column at fixed row_idx, vary along row_dim
       If none of 'names' found/matching dims, return ones."""
    for nm in names:
        if nm in ds:
            da = ds[nm]
            if 'time' in da.dims:
                da = da.mean('time')
            dims = da.dims
            if (row_dim in dims) and (col_dim in dims):
                if prefer_row:
                    line = da.isel({row_dim: row_idx}).transpose(..., col_dim)
                    arr  = np.asarray(line.values, dtype=float)
                else:
                    line = da.isel({col_dim: row_idx}).transpose(..., row_dim)
                    arr  = np.asarray(line.values, dtype=float)
                if np.all(~np.isfinite(arr)):
                    continue
                # replace NaNs by the finite median to keep spacing usable
                finite = np.isfinite(arr)
                if finite.any():
                    med = np.nanmedian(arr)
                    arr[~finite] = med
                else:
                    arr[:] = 1.0
                return arr
    # fallback: uniform spacing = 1 cell
    return None  # signal "uniform cells"

def positions_from_spacing(spacing, n):
    if spacing is None:
        return np.arange(n, dtype=float)  # 0,1,2,...
    # spacing[k] ≈ distance between k and k+1; build 0-based positions
    s = np.asarray(spacing, dtype=float)
    s = s[:n-1] if s.size >= (n-1) else np.pad(s, (0, n-1-s.size), constant_values=float(np.nanmedian(s)))
    pos = np.zeros(n, dtype=float)
    pos[1:] = np.cumsum(s)
    return pos

def grad1d(y, x):
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    g = np.full_like(y, np.nan)
    n = y.size
    # finite mask
    fin = np.isfinite(y)
    if fin.sum() < 3:
        return g  # not enough points
    # interior points: centered differences where neighbors finite
    for i in range(n):
        if i == 0:
            # forward difference if next is finite
            j = 1
            if fin[i] and fin[j]:
                dx = x[j] - x[i]
                if dx != 0: g[i] = (y[j] - y[i]) / dx
        elif i == n-1:
            j = n-2
            if fin[i] and fin[j]:
                dx = x[i] - x[j]
                if dx != 0: g[i] = (y[i] - y[j]) / dx
        else:
            i0, i1 = i-1, i+1
            if fin[i0] and fin[i1]:
                dx = x[i1] - x[i0]
                if dx != 0: g[i] = (y[i1] - y[i0]) / dx
            elif fin[i] and fin[i1]:
                dx = x[i1] - x[i]
                if dx != 0: g[i] = (y[i1] - y[i]) / dx
            elif fin[i] and fin[i0]:
                dx = x[i] - x[i0]
                if dx != 0: g[i] = (y[i] - y[i0]) / dx
    return g

def time_derivative_face(da_inst):
    """Centered finite-difference d(da)/dt along time (seconds).
       At endpoints use one-sided; returns same shape as input."""
    tvals = da_inst['time'].values
    # seconds since start
    tsec = hours_since_start(tvals) * 3600.0
    # build per-interval dt
    dt = np.diff(tsec)  # length Nt-1
    # allocate output
    darr = da_inst.values
    deriv = np.empty_like(darr)
    # endpoints: one-sided
    deriv[0,...]  = (darr[1,...] - darr[0,...]) / (dt[0] if dt[0] != 0 else 1.0)
    deriv[-1,...] = (darr[-1,...] - darr[-2,...]) / (dt[-1] if dt[-1] != 0 else 1.0)
    # interior: centered
    for k in range(1, len(tsec)-1):
        dtc = tsec[k+1] - tsec[k-1]
        deriv[k,...] = (darr[k+1,...] - darr[k-1,...]) / (dtc if dtc != 0 else 1.0)
    return xr.DataArray(deriv, dims=da_inst.dims, coords=da_inst.coords, name=f"d{da_inst.name}_dt")

# -------------------
# LOAD DATA
# -------------------
daily = load_daily()
inst  = load_inst()

# face vars (instantaneous)
uE_name = pick_var(inst, 'uvelE', ['uvelE_1','uvelE_0001','uvelE_01'])
vN_name = pick_var(inst, 'vvelN', ['vvelN_1','vvelN_0001','vvelN_01'])
uE_inst = inst[uE_name]        # (time, nj, ni_E)
vN_inst = inst[vN_name]        # (time, nj_N, ni)

# daily mean face stresses + T-grid speeds
kuxe_name = pick_var(daily, 'KuxE')
kuyn_name = pick_var(daily, 'KuyN')
uT_name   = pick_var(daily, 'uvel')
vT_name   = pick_var(daily, 'vvel')

KuxE = daily[kuxe_name]        # (time, nj, ni)
KuyN = daily[kuyn_name]
uT   = daily[uT_name]
vT   = daily[vT_name]

# edge masks
E_mask, N_mask = build_edge_masks(daily)  # boolean (nj,ni)

# get dim names for pretty plotting
nj_name_E, ni_name_E = dims_ij(uE_inst.isel(time=0))
nj_name_N, ni_name_N = dims_ij(vN_inst.isel(time=0))

# -------------------
# STATS
# -------------------
# choose last instantaneous output as "representative snapshot"
uE_last = uE_inst.isel(time=-1)
vN_last = vN_inst.isel(time=-1)

# edge vs interior (instantaneous faces)
uE_edge_mean  = masked_mean_abs(uE_last, E_mask)
uE_int_mean   = masked_mean_abs(uE_last, ~E_mask)
vN_edge_mean  = masked_mean_abs(vN_last, N_mask)
vN_int_mean   = masked_mean_abs(vN_last, ~N_mask)

# daily means
KuxE_edge = float(np.abs(KuxE).where(E_mask).mean().values)
KuxE_int  = float(np.abs(KuxE).where(~E_mask).mean().values)
KuyN_edge = float(np.abs(KuyN).where(N_mask).mean().values)
KuyN_int  = float(np.abs(KuyN).where(~N_mask).mean().values)

# time series of spatial mean |uE| over edge vs interior
uE_edge_ts = np.abs(uE_inst).where(E_mask).mean(dim=(nj_name_E, ni_name_E))
uE_int_ts  = np.abs(uE_inst).where(~E_mask).mean(dim=(nj_name_E, ni_name_E))
vN_edge_ts = np.abs(vN_inst).where(N_mask).mean(dim=(nj_name_N, ni_name_N))
vN_int_ts  = np.abs(vN_inst).where(~N_mask).mean(dim=(nj_name_N, ni_name_N))

# T-grid “edge” = anywhere CDP is active on adjacent faces
T_edge_mask = (E_mask | N_mask)

# Pick vars (prefer T-grid versions; fall back to E/N if needed)
aice_name = pick_var(daily, 'aice', [])
hi_name   = pick_var(daily, 'hi',   [])

try:
    strx_name = pick_var(daily, 'strintx', ['strintxE'])
except KeyError:
    strx_name = 'strintxE'  # if only E exists
try:
    stry_name = pick_var(daily, 'strinty', ['strintyN'])
except KeyError:
    stry_name = 'strintyN'  # if only N exists

AICE = daily[aice_name]
HI   = daily[hi_name]
STRX = daily[strx_name]
STRY = daily[stry_name]

# Compute edge vs interior means
aice_edge = mean_over(T_edge_mask, AICE, use_abs=False)
aice_int  = mean_over(~T_edge_mask, AICE, use_abs=False)
hi_edge   = mean_over(T_edge_mask, HI,   use_abs=False)
hi_int    = mean_over(~T_edge_mask, HI,  use_abs=False)

# stresses: use absolute value (magnitude)
# If STRX/STRY are on E/N faces, the broadcast will align with those shapes automatically.
strx_edge = mean_over(E_mask if 'E' in strx_name[-1:] else T_edge_mask, STRX, use_abs=True)
strx_int  = mean_over(~(E_mask if 'E' in strx_name[-1:] else T_edge_mask), STRX, use_abs=True)
stry_edge = mean_over(N_mask if 'N' in stry_name[-1:] else T_edge_mask, STRY, use_abs=True)
stry_int  = mean_over(~(N_mask if 'N' in stry_name[-1:] else T_edge_mask), STRY, use_abs=True)

# time derivatives
duE_dt = time_derivative_face(uE_inst)
dvN_dt = time_derivative_face(vN_inst)

# edge vs interior mean |d/dt|
uE_dt_edge = float(np.abs(duE_dt).where(E_mask).mean(dim=('time', nj_name_E, ni_name_E)).values)
uE_dt_int  = float(np.abs(duE_dt).where(~E_mask).mean(dim=('time', nj_name_E, ni_name_E)).values)
vN_dt_edge = float(np.abs(dvN_dt).where(N_mask).mean(dim=('time', nj_name_N, ni_name_N)).values)
vN_dt_int  = float(np.abs(dvN_dt).where(~N_mask).mean(dim=('time', nj_name_N, ni_name_N)).values)

# ---------- center-line spatial gradients on daily T-grid means ----------
# 1) Build daily-mean face velocities from instantaneous files
uE_mean = uE_inst.mean('time')    # (nj_E, ni_E)
vN_mean = vN_inst.mean('time')    # (nj_N, ni_N)

# 2) Center-line indices on each face grid
nyE = uE_mean.sizes[nj_name_E]; nxE = uE_mean.sizes[ni_name_E]
nyN = vN_mean.sizes[nj_name_N]; nxN = vN_mean.sizes[ni_name_N]
jmid_E = nyE // 2
imid_N = nxN // 2

# 3) Get spacings along the center lines
dx_line = get_spacing_line(daily, names=('dxE','hte','HTE','dxu'),
                           row_dim=nj_name_E, col_dim=ni_name_E,
                           row_idx=jmid_E, prefer_row=True)   # along i on E
dy_line = get_spacing_line(daily, names=('dyN','htn','HTN','dyv'),
                           row_dim=nj_name_N, col_dim=ni_name_N,
                           row_idx=imid_N, prefer_row=False)  # along j on N

xE = positions_from_spacing(dx_line, nxE)
yN = positions_from_spacing(dy_line, nyN)
x_units = 'm' if dx_line is not None else 'cell'
y_units = 'm' if dy_line is not None else 'cell'

# 4) Extract center-line profiles on faces
uE_center = uE_mean.isel({nj_name_E: jmid_E}).values  # shape (nxE,)
vN_center = vN_mean.isel({ni_name_N: imid_N}).values  # shape (nyN,)

duE_dx = grad1d(uE_center, xE)   # (nxE,)
dvN_dy = grad1d(vN_center, yN)   # (nyN,)

# 5) Stats from middle → edges (ignore NaNs)
# define slices
u_right  = duE_dx[nxE//2 : ]
u_left   = duE_dx[: nxE//2 + 1]
v_top    = dvN_dy[nyN//2 : ]
v_bottom = dvN_dy[: nyN//2 + 1]

du_dx_mean_R = float(np.nanmean(np.abs(u_right)))
du_dx_mean_L = float(np.nanmean(np.abs(u_left )))
dv_dy_mean_T = float(np.nanmean(np.abs(v_top  )))
dv_dy_mean_B = float(np.nanmean(np.abs(v_bottom)))

print('--- CDP diagnostics (', SIM_LABEL, ') ---')
print(f'Edge vs Interior (instantaneous, last time):')
print(f'  <|uE|> edge = {uE_edge_mean:9.4e}   interior = {uE_int_mean:9.4e}   ratio edge/int = {uE_edge_mean/(uE_int_mean+EPS):.3f}')
print(f'  <|vN|> edge = {vN_edge_mean:9.4e}   interior = {vN_int_mean:9.4e}   ratio edge/int = {vN_edge_mean/(vN_int_mean+EPS):.3f}')
print('Time-derivative stats (instantaneous faces):')
print(f"  <|duE/dt|> edge = {uE_dt_edge:9.4e}   interior = {uE_dt_int:9.4e}")
print(f"  <|dvN/dt|> edge = {vN_dt_edge:9.4e}   interior = {vN_dt_int:9.4e}")
print('Center-line spatial gradients on faces:')
print(f"  mean |∂u/∂x| (E face, mid→east) = {du_dx_mean_R:9.4e}   (mid→west) = {du_dx_mean_L:9.4e}   [{x_units}⁻¹]")
print(f"  mean |∂v/∂y| (N face, mid→north)= {dv_dy_mean_T:9.4e}   (mid→south) = {dv_dy_mean_B:9.4e}   [{y_units}⁻¹]")
print(f'Daily-mean stresses:')
print(f'  <|KuxE|> edge = {KuxE_edge:9.4e}  interior = {KuxE_int:9.4e}')
print(f'  <|KuyN|> edge = {KuyN_edge:9.4e}  interior = {KuyN_int:9.4e}')
print('Extra daily stats (T-grid unless noted):')
print(f"  <aice>  edge = {aice_edge:9.4e}   interior = {aice_int:9.4e}   Δ(edge-int) = {aice_edge-aice_int: .2e}")
print(f"  <hi>    edge = {hi_edge:9.4e}   interior = {hi_int:9.4e}   Δ(edge-int) = {hi_edge-hi_int: .2e}")
print(f"  <|strintx|> edge = {strx_edge:9.4e}   interior = {strx_int:9.4e}")
print(f"  <|strinty|> edge = {stry_edge:9.4e}   interior = {stry_int:9.4e}")

# =======================
# Assemble ALL stats into a CSV
# =======================
# From earlier (instantaneous faces & daily stresses)
push('uE_inst_last', 'E', uE_edge_mean, uE_int_mean, use_abs=True, note='instantaneous |uE| at last time')
push('vN_inst_last', 'N', vN_edge_mean, vN_int_mean, use_abs=True, note='instantaneous |vN| at last time')
push('KuxE_daily',   'E', KuxE_edge,    KuxE_int,    use_abs=True, note='daily mean |KuxE|')
push('KuyN_daily',   'N', KuyN_edge,    KuyN_int,    use_abs=True, note='daily mean |KuyN|')
# New T-grid vars
push('aice_daily',   'T', aice_edge, aice_int, use_abs=False)
push('hi_daily',     'T', hi_edge,   hi_int,   use_abs=False)
# Internal stress divergence components (may be T or E/N depending on output)
grid_x = 'E' if strx_name.lower().endswith('e') else 'T'
grid_y = 'N' if stry_name.lower().endswith('n') else 'T'
push('strintx_daily', grid_x, strx_edge, strx_int, use_abs=True)
push('strinty_daily', grid_y, stry_edge, stry_int, use_abs=True)
push('duE_dt_inst', 'E', uE_dt_edge, uE_dt_int, use_abs=True, note='mean |duE/dt| over time')
push('dvN_dt_inst', 'N', vN_dt_edge, vN_dt_int, use_abs=True, note='mean |dvN/dt| over time')
# store center-line gradient summary (as pseudo “edge” vs “interior” columns)
push('du_dx_mid→east', 'T', du_dx_mean_R, 0.0, use_abs=True, note=f'center-line; units per {x_units}')
push('du_dx_mid→west', 'T', du_dx_mean_L, 0.0, use_abs=True, note=f'center-line; units per {x_units}')
push('dv_dy_mid→north','T', dv_dy_mean_T, 0.0, use_abs=True, note=f'center-line; units per {y_units}')
push('dv_dy_mid→south','T', dv_dy_mean_B, 0.0, use_abs=True, note=f'center-line; units per {y_units}')
df = pd.DataFrame(rows, columns=['sim','var','grid','abs','edge_mean','interior_mean','ratio','note'])
csv_path = os.path.join(outdir, f'cdp_stats_{SIM_LABEL}.csv')
df.to_csv(csv_path, index=False)
print(f"\nCSV written: {csv_path}")

# -------------------
# FIGURES
# -------------------
# Fig 1: instantaneous uE snapshot + edge mask outline
plt.figure(figsize=(8,6))
im = plt.imshow(uE_last.values, origin='lower', aspect='equal')
plt.colorbar(im, label='uE (m/s)')
# overlay mask (thin lines): draw mask as contour if available
try:
    plt.contour(E_mask.values.astype(float), levels=[0.5], colors='w', linewidths=1.0)
except Exception:
    pass
plt.title(f'Instantaneous uE (last output) — {SIM_LABEL}')
plt.xlabel(ni_name_E); plt.ylabel(nj_name_E)
f1 = os.path.join(outdir, f'fig1_uE_inst_{SIM_LABEL}.png'); plt.savefig(f1, dpi=150, bbox_inches='tight')#; plt.show()

# Fig 2: daily mean |KuxE|
plt.figure(figsize=(8,6))
im = plt.imshow(np.abs(KuxE.mean('time').values), origin='lower', aspect='equal')
plt.colorbar(im, label='|KuxE| (N m$^{-2}$)')
plt.contour(E_mask.values.astype(float), levels=[0.5], colors='k', linewidths=0.8)
plt.title(f'Daily mean |KuxE| — {SIM_LABEL}')
plt.xlabel(ni_name_E); plt.ylabel(nj_name_E)
f2 = os.path.join(outdir, f'fig2_absKuxE_daily_{SIM_LABEL}.png'); plt.savefig(f2, dpi=150, bbox_inches='tight')#; plt.show()

# Fig 3: daily mean |KuyN|
plt.figure(figsize=(8,6))
im = plt.imshow(np.abs(KuyN.mean('time').values), origin='lower', aspect='equal')
plt.colorbar(im, label='|KuyN| (N m$^{-2}$)')
plt.contour(N_mask.values.astype(float), levels=[0.5], colors='k', linewidths=0.8)
plt.title(f'Daily mean |KuyN| — {SIM_LABEL}')
nj_name_N, ni_name_N = dims_ij(KuyN.isel(time=0))
plt.xlabel(ni_name_N); plt.ylabel(nj_name_N)
f3 = os.path.join(outdir, f'fig3_absKuyN_daily_{SIM_LABEL}.png'); plt.savefig(f3, dpi=150, bbox_inches='tight')#; plt.show()

# Fig 4: time series of |uE| and |vN| (edge vs interior)
# Convert CFTime to numeric hours since start (works for 360_day)

t_hours_u = hours_since_start(uE_edge_ts['time'].values)
t_hours_v = hours_since_start(vN_edge_ts['time'].values)
fig, ax = plt.subplots(2,1, figsize=(10,6), sharex=True)
ax[0].plot(t_hours_u, uE_edge_ts.values, label='edge')
ax[0].plot(t_hours_u, uE_int_ts.values , label='interior', linestyle='--')
ax[0].set_ylabel('<|uE|> (m/s)'); ax[0].legend(); ax[0].set_title('Instantaneous |uE| (spatial mean)')
ax[1].plot(t_hours_v, vN_edge_ts.values, label='edge')
ax[1].plot(t_hours_v, vN_int_ts.values , label='interior', linestyle='--')
ax[1].set_ylabel('<|vN|> (m/s)'); ax[1].legend(); ax[1].set_title('Instantaneous |vN| (spatial mean)')
ax[1].set_xlabel('hours since start')
f4 = os.path.join(outdir, f'fig4_timeseries_{SIM_LABEL}.png')
plt.savefig(f4, dpi=150, bbox_inches='tight')
plt.close(fig)

# Fig U (E face): u profile & ∂u/∂x
fig, ax = plt.subplots(1,2, figsize=(12,4))
ax[0].plot(xE, uE_center)
ax[0].set_title('u (E face) center-line (j=mid)')
ax[0].set_xlabel(f'x [{x_units}]'); ax[0].set_ylabel('u (m s$^{-1}$)')
ax[1].plot(xE, duE_dx)
ax[1].set_title('∂u/∂x (E face, center-line)')
ax[1].set_xlabel(f'x [{x_units}]'); ax[1].set_ylabel(f'∂u/∂x (m s$^{{-1}}$ per {x_units})')
f5 = os.path.join(outdir, f'fig5E_centerline_u_du_dx_{SIM_LABEL}.png')
plt.savefig(f5, dpi=150, bbox_inches='tight'); plt.close(fig)

# Fig V (N face): v profile & ∂v/∂y
fig, ax = plt.subplots(1,2, figsize=(12,4))
ax[0].plot(yN, vN_center)
ax[0].set_title('v (N face) center-line (i=mid)')
ax[0].set_xlabel(f'y [{y_units}]'); ax[0].set_ylabel('v (m s$^{-1}$)')
ax[1].plot(yN, dvN_dy)
ax[1].set_title('∂v/∂y (N face, center-line)')
ax[1].set_xlabel(f'y [{y_units}]'); ax[1].set_ylabel(f'∂v/∂y (m s$^{{-1}}$ per {y_units})')
f6 = os.path.join(outdir, f'fig6N_centerline_v_dv_dy_{SIM_LABEL}.png')
plt.savefig(f6, dpi=150, bbox_inches='tight'); plt.close(fig)

print("Saved:")
for p in (f1,f2,f3,f4,f5,f6):
    print("  ", p)
