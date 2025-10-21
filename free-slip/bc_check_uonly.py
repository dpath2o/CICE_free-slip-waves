#!/usr/bin/env python3
import argparse, glob, os
import numpy as np
import xarray as xr

def latest_inst(path):
    files = sorted(glob.glob(os.path.join(path, "iceh_inst.*.nc")))
    if not files: raise FileNotFoundError(f"No iceh_inst* in {path}")
    return files[-1]

def get_ds_file(path, fname=None):
    f = fname if fname else latest_inst(path)
    ds = xr.open_dataset(f, decode_cf=True, mask_and_scale=True)
    return ds, f

def arr(ds, name, sel=None):
    """time slice -> np.array with NaNs for missing/fill"""
    if name not in ds: return None
    a = ds[name]
    if "time" in a.dims:
        a = a.isel(time=0)
    if sel is not None:
        a = a.sel(**sel)
    return a.to_masked_array().filled(np.nan).astype(np.float64)

def safe_bool_from_mask(a, thresh=0.5):
    """tmask/umask -> boolean; NaN treated as False (land)"""
    return np.where(np.isfinite(a), a > thresh, False)

def p99_abs(a, m):
    v = np.abs(a[np.logical_and(m, np.isfinite(a))])
    if v.size == 0: return np.nan, 0
    return np.nanpercentile(v, 99.0), v.size

def print_val(label, val, nvalid, scale=None):
    if np.isnan(val):
        print(f"  {label}: n/a (no finite coastal ocean samples)")
    else:
        if scale is not None and not np.isnan(scale) and scale > 0:
            print(f"  {label}: {val: .3e} (rel={val/scale: .3e}, n={nvalid})")
        else:
            print(f"  {label}: {val: .3e} (n={nvalid})")

def build_coast_masks_from_tmask(tmask):
    """
    tmask (nj,ni) -> coastE (nj,ni-1), coastN (nj-1,ni)
    A face is coastal if the two adjacent T cells differ in ocean/land.
    """
    tocean = safe_bool_from_mask(tmask)
    coastE = np.logical_xor(tocean[:, :-1], tocean[:, 1:])
    coastN = np.logical_xor(tocean[:-1, :], tocean[1:, :])
    return coastE, coastN, tocean

def map_faces_to_U(coastE, coastN, umask):
    """
    Map E/N coastal faces to nearby U corners; keep only U with umask==1.
    For an E face between T(i,j) and T(i+1,j), adjacent U are (i+1,j) and (i+1,j+1).
    For an N face between T(i,j) and T(i,j+1), adjacent U are (i,j+1) and (i+1,j+1).
    """
    nj, ni = umask.shape
    UE = np.zeros((nj, ni), dtype=bool)
    UN = np.zeros((nj, ni), dtype=bool)

    # E-face -> U(i+1,j) and U(i+1,j+1)
    if coastE.size:
        UE[:, 1:] |= coastE            # (j, i+1)
        UE[1:, 1:] |= coastE[:-1, :]   # (j+1, i+1)

    # N-face -> U(i,j+1) and U(i+1,j+1)
    if coastN.size:
        UN[1:, :] |= coastN            # (j+1, i)
        UN[1:, 1:] |= coastN[:, :-1]   # (j+1, i+1)

    Uocean = safe_bool_from_mask(umask)
    UE &= Uocean
    UN &= Uocean
    return UE, UN

def interior_scale(vE_tan, uN_tan, UE, UN, umask):
    Uocean = safe_bool_from_mask(umask)
    intE = np.logical_and(Uocean, ~UE)
    intN = np.logical_and(Uocean, ~UN)
    sE, nE = p99_abs(vE_tan, intE)
    sN, nN = p99_abs(uN_tan, intN)
    return sE, sN

def one_sided_derivs(u, v, dxu, dyu):
    """
    Compute forward differences (one-sided) with NaN safety.
    dv/dx at U: (v[j,i] - v[j,i-1]) / dxu[j,i], placed at i>=1
    du/dy at U: (u[j,i] - u[j-1,i]) / dyu[j,i], placed at j>=1
    If either operand or metric is NaN/zero, result is NaN.
    """
    nj, ni = u.shape
    dv_dx = np.full_like(v, np.nan)
    du_dy = np.full_like(u, np.nan)

    # dv/dx
    denom_x = dxu.copy() if dxu is not None else np.ones_like(v)
    denom_x = np.where((np.isfinite(denom_x)) & (denom_x != 0.0), denom_x, np.nan)
    dv = v[:, 1:] - v[:, :-1]
    dv_dx[:, 1:] = dv / denom_x[:, 1:]

    # du/dy
    denom_y = dyu.copy() if dyu is not None else np.ones_like(u)
    denom_y = np.where((np.isfinite(denom_y)) & (denom_y != 0.0), denom_y, np.nan)
    du = u[1:, :] - u[:-1, :]
    du_dy[1:, :] = du / denom_y[1:, :]

    return dv_dx, du_dy

def check_one(ds, label, fname=None):
    ds, used = get_ds_file(ds, fname)
    # Required fields
    tmask = arr(ds, "tmask")
    umask = arr(ds, "umask")
    u = arr(ds, "uvel_1")
    v = arr(ds, "vvel_1")
    dxu = arr(ds, "dxu")
    dyu = arr(ds, "dyu")

    # Basic existence checks
    for n, A in [("tmask", tmask), ("umask", umask), ("uvel_1", u), ("vvel_1", v)]:
        if A is None:
            raise RuntimeError(f"Missing variable {n} in {used}")

    coastE, coastN, _ = build_coast_masks_from_tmask(tmask)
    UE, UN = map_faces_to_U(coastE, coastN, umask)

    # Tangential on coasts (at U): E wall -> v ; N wall -> u
    vE_tan = v
    uN_tan = u

    # Interior scales away from coast
    sE, sN = interior_scale(vE_tan, uN_tan, UE, UN, umask)

    # Normals at coasts (at U): E wall -> u ; N wall -> v
    p99_unE, n_unE = p99_abs(u, UE)
    p99_vnN, n_vnN = p99_abs(v, UN)

    # Tangentials at coasts
    p99_vtanE, n_vtanE = p99_abs(v, UE)
    p99_utanN, n_utanN = p99_abs(u, UN)

    # One-sided normal derivatives of tangential
    dv_dx, du_dy = one_sided_derivs(u, v, dxu, dyu)
    p99_dvdx_E, n_dvdx_E = p99_abs(dv_dx, UE)
    p99_dudy_N, n_dudy_N = p99_abs(du_dy, UN)

    # Summaries
    print(f"\n=== {label} ({os.path.basename(used)}) ===")
    print(f"Coast faces: E={int(coastE.sum())}, N={int(coastN.sum())}; coastal U: UE={int(UE.sum())}, UN={int(UN.sum())}")
    print(f"Interior scales (99th pct away from coast): V={sE: .3e} (E tangential), U={sN: .3e} (N tangential)")

    print("Normals ~ 0 (both BCs):")
    print_val("p99 |u| on E-coast-U", p99_unE, n_unE, sN)  # use U scale
    print_val("p99 |v| on N-coast-U", p99_vnN, n_vnN, sE)  # use V scale

    print("Tangentials on coast (no-slip -> ~0; free-slip unconstrained):")
    print_val("p99 |v| on E-coast-U", p99_vtanE, n_vtanE, sE)
    print_val("p99 |u| on N-coast-U", p99_utanN, n_utanN, sN)

    print("Normal derivative of tangential (free-slip -> ~0):")
    print_val("max |∂v/∂x| on E-coast-U", np.nanmax(np.abs(dv_dx[UE])) if UE.any() else np.nan, n_dvdx_E)
    print_val("max |∂u/∂y| on N-coast-U", np.nanmax(np.abs(du_dy[UN])) if UN.any() else np.nan, n_dudy_N)

    # Verdicts (graceful if n/a)
    def ok(val, scale, tol_rel=5e-3, tol_abs=5e-6):
        if np.isnan(val): return True  # nothing to test -> don't fail hard
        # interpret in "m/s" for velocities, and 1/m for derivatives; tol_abs is small
        thr = max(tol_abs, tol_rel * (scale if (scale is not None and not np.isnan(scale)) else 1.0))
        return val <= thr

    exp = label.lower()
    is_noslip = ("noslip" in exp)
    is_free   = ("free" in exp) or ("noratio" in exp)

    imp_ok = ok(p99_unE, sN) and ok(p99_vnN, sE)
    if is_noslip:
        tan_ok = ok(p99_vtanE, sE) and ok(p99_utanN, sN)
        print(f"Verdict: impermeability={'PASS' if imp_ok else 'FAIL'}, no-slip tangential={'PASS' if tan_ok else 'FAIL'}")
    elif is_free:
        d_ok = ok(np.nanmax(np.abs(dv_dx[UE])) if UE.any() else np.nan, 1.0) and \
               ok(np.nanmax(np.abs(du_dy[UN])) if UN.any() else np.nan, 1.0)
        print(f"Verdict: impermeability={'PASS' if imp_ok else 'FAIL'}, free-slip d/dn tangential={'PASS' if d_ok else 'FAIL'}")
    else:
        print(f"Verdict: impermeability={'PASS' if imp_ok else 'FAIL'}")

def main():
    ap = argparse.ArgumentParser(description="BC checks using only U-point velocities")
    ap.add_argument("--exp", action="append", required=True,
                    help="label:/path/to/history (label used to infer BC type: contains 'noslip' or 'free'/'noratio')")
    ap.add_argument("--file", default=None, help="specific iceh_inst.*.nc filename (basename only); default: latest")
    args = ap.parse_args()

    for item in args.exp:
        if ":" not in item:
            raise SystemExit("--exp must be like name:/path")
        label, path = item.split(":", 1)
        fname = None
        if args.file:
            fname = os.path.join(path, args.file)
        check_one(path, label, fname)

if __name__ == "__main__":
    main()
