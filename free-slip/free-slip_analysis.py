#!/usr/bin/env python3
import numpy as np
import xarray as xr
import os, glob

def pick(ds, *names):
    for n in names:
        if n in ds.variables:
            return ds[n]
    raise KeyError(f"None of {names} found")

def robust_stats(arr, mask, label=""):
    a = np.where(mask, arr, np.nan)
    a = np.abs(a[np.isfinite(a)])
    if a.size == 0:
        return dict(count=0, max=0.0, mean=0.0, p95=0.0, note=f"{label}: no points")
    return dict(count=a.size, max=float(a.max()), mean=float(a.mean()),
                p95=float(np.percentile(a, 95)), note="")

def load_one(path):
    files = sorted(glob.glob(os.path.join(path, "iceh_inst.*.nc")))
    if not files:
        raise FileNotFoundError(f"No inst files in {path}")
    return files, xr.open_dataset(files[-1], decode_times=False)

def build_coast_masks(ds):
    """
    E/N coastal faces:
    1) Prefer F2E_1 / F2N_1 > 0 (coastal form factors in the file)
    2) Fallback: domain perimeter (i==0 or i==ni-1 for E; j==0 or j==nj-1 for N),
       filtered by emask/nmask > 0.5
    U-corner coast = any touching coastal E/N faces.
    """
    nj, ni = pick(ds, "tmask").shape  # (nj, ni)

    # Try coastal form factors
    have_F2E = "F2E_1" in ds.variables
    have_F2N = "F2N_1" in ds.variables
    if have_F2E:
        coastE = np.asarray(ds["F2E_1"][0, :, :]) > 0
    else:
        coastE = np.zeros((nj, ni), dtype=bool)
    if have_F2N:
        coastN = np.asarray(ds["F2N_1"][0, :, :]) > 0
    else:
        coastN = np.zeros((nj, ni), dtype=bool)

    # Fallback to perimeter if needed (or if the form factors yield 0 points)
    if not coastE.any() or not coastN.any():
        em = pick(ds, "emask").values > 0.5
        nm = pick(ds, "nmask").values > 0.5
        if not coastE.any():
            coastE = em & (np.equal.outer(np.arange(nj), np.arange(nj)) * 0 + 1).astype(bool)  # dummy to get shape
            coastE[:] = False
            # West/East “wall” columns (faces live on same (nj,ni) grid in outputs)
            westE = np.zeros_like(em, dtype=bool);  westE[:, 0]     = True
            eastE = np.zeros_like(em, dtype=bool);  eastE[:, ni-1]  = True
            coastE = (westE | eastE) & em
        if not coastN.any():
            coastN = nm & (np.equal.outer(np.arange(nj), np.arange(nj)) * 0 + 1).astype(bool)
            coastN[:] = False
            southN = np.zeros_like(nm, dtype=bool); southN[0, :]     = True
            northN = np.zeros_like(nm, dtype=bool); northN[nj-1, :]  = True
            coastN = (southN | northN) & nm

    # U-corner is coastal if any of its four adjacent faces is coastal
    coastU = np.zeros((nj, ni), dtype=bool)
    cornerU = np.zeros((nj, ni), dtype=bool)
    for j in range(1, nj-1):
        for i in range(1, ni-1):
            w = coastE[j, i-1]
            e = coastE[j, i  ]
            s = coastN[j-1, i]
            n = coastN[j,   i]
            c = w or e or s or n
            coastU[j, i]  = c
            cornerU[j, i] = (w and s) or (w and n) or (e and s) or (e and n)

    return coastE, coastN, coastU, cornerU

def normal_derivative_on_E(ds, coastE, vE):
    """One-sided index-space normal derivative of v (tangential on E) at E-coast."""
    tmB = pick(ds, "tmask").values > 0.5
    nj, ni = tmB.shape
    deriv = np.full_like(vE, np.nan)
    js, is_ = np.where(coastE)
    for j, i in zip(js, is_):
        # If interior to the +x side, use (i+1) – i, else i – (i-1)
        if i+1 < ni and tmB[j, min(i+1, ni-1)]:
            ni2 = i+1
            deriv[j, i] = vE[j, ni2] - vE[j, i]
        elif i-1 >= 0 and tmB[j, max(i-1, 0)]:
            ni2 = i-1
            deriv[j, i] = vE[j, i] - vE[j, ni2]
    return deriv

def normal_derivative_on_N(ds, coastN, uN):
    """One-sided index-space normal derivative of u (tangential on N) at N-coast."""
    tmB = pick(ds, "tmask").values > 0.5
    nj, ni = tmB.shape
    deriv = np.full_like(uN, np.nan)
    js, is_ = np.where(coastN)
    for j, i in zip(js, is_):
        if j+1 < nj and tmB[min(j+1, nj-1), i]:
            nj2 = j+1
            deriv[j, i] = uN[nj2, i] - uN[j, i]
        elif j-1 >= 0 and tmB[max(j-1, 0), i]:
            nj2 = j-1
            deriv[j, i] = uN[j, i] - uN[nj2, i]
    return deriv

def interior_scale(v, mask):
    a = np.abs(np.where(mask, v, np.nan))
    a = a[np.isfinite(a)]
    return float(np.percentile(a, 99)) if a.size else 1.0  # avoid div-by-zero

def analyze_snapshot(path, label):
    files, ds = load_one(path)
    coastE, coastN, coastU, _ = build_coast_masks(ds)

    # velocities (single time sample)
    uE = pick(ds, "uvelE_1")[0].values
    vE = pick(ds, "vvelE_1")[0].values
    uN = pick(ds, "uvelN_1")[0].values
    vN = pick(ds, "vvelN_1")[0].values

    em = pick(ds, "emask").values > 0.5
    nm = pick(ds, "nmask").values > 0.5
    interior_E = em & ~coastE
    interior_N = nm & ~coastN

    Vref_E = interior_scale(vE, interior_E)
    Uref_N = interior_scale(uN, interior_N)

    # normals (both BCs should be ~0)
    stat_un_E = robust_stats(uE, coastE, "u_n(E)")
    stat_vn_N = robust_stats(vN, coastN, "v_n(N)")

    # tangentials (no-slip ~0; free-slip can be non-zero)
    stat_vt_E = robust_stats(vE, coastE, "v_tan(E)")
    stat_ut_N = robust_stats(uN, coastN, "u_tan(N)")

    # normal derivatives of tangential (free-slip ~0)
    dvE = normal_derivative_on_E(ds, coastE, vE)
    duN = normal_derivative_on_N(ds, coastN, uN)
    stat_dv_dn_E = robust_stats(dvE, np.isfinite(dvE), "∂v/∂n(E)")
    stat_du_dn_N = robust_stats(duN, np.isfinite(duN), "∂u/∂n(N)")

    # shear at coastal U (diagnostic)
    if "shearU_1" in ds:
        shearU = pick(ds, "shearU_1")[0].values
        stat_shearU = robust_stats(shearU, coastU, "shearU(Ucoast)")
    else:
        stat_shearU = dict(count=0, max=0.0, mean=0.0, p95=0.0, note="no shearU_1")

    print(f"\n=== {label} ({os.path.basename(files[-1])}) ===")
    print(f"Coast E faces: {stat_un_E['count']}, Coast N faces: {stat_vn_N['count']}")
    if stat_un_E['count'] == 0 or stat_vn_N['count'] == 0:
        print("  [info] No coastal faces were tagged; used perimeter fallback. If you expected tagged coast,")
        print("         check F2E_1/F2N_1 fields or confirm this is a closed box with walls not stored as land.")
    print(f"Interior scales (99th pct): V_E={Vref_E:.3e}, U_N={Uref_N:.3e}")

    print("Normals ~ 0 (both BCs):")
    print(f"  max |u_n| on E coast   : {stat_un_E['max']:.3e}")
    print(f"  max |v_n| on N coast   : {stat_vn_N['max']:.3e}")

    print("Tangentials (no-slip ~0; free-slip unconstrained):")
    rel_vt = stat_vt_E['max']/Vref_E if Vref_E>0 else 0.0
    rel_ut = stat_ut_N['max']/Uref_N if Uref_N>0 else 0.0
    print(f"  max |v_tan| on E coast : {stat_vt_E['max']:.3e}  (rel={rel_vt:.3e})")
    print(f"  max |u_tan| on N coast : {stat_ut_N['max']:.3e}  (rel={rel_ut:.3e})")

    print("Normal derivative of tangential (free-slip ~0):")
    print(f"  max |dv/dn| on E coast : {stat_dv_dn_E['max']:.3e} (index units)")
    print(f"  max |du/dn| on N coast : {stat_du_dn_N['max']:.3e} (index units)")

    print("Shear on coastal U (diagnostic):")
    print(f"  max |shearU| at coastal U : {stat_shearU['max']:.3e}")


def main():
    free_dir   = "/g/data/gv90/da1339/cice-dirs/runs/BCtest1_free/history/"
    noslip_dir = "/g/data/gv90/da1339/cice-dirs/runs/BCtest1_noslip/history/"
    analyze_snapshot(free_dir,   "FREE-SLIP")
    analyze_snapshot(noslip_dir, "NO-SLIP")

if __name__ == "__main__":
    main()