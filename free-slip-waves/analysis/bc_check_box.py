
#!/usr/bin/env python3
import argparse, re, glob, os, sys
from datetime import datetime, date
from typing import List, Tuple, Optional, Dict
import numpy as np
import xarray as xr

FN_RE = re.compile(r".*iceh_inst\.(\d{4})-(\d{2})-(\d{2})-(\d{5})\.nc$")

def parse_datestr(s: str) -> date:
    return datetime.strptime(s, "%Y-%m-%d").date()

def file_date(f: str) -> Optional[date]:
    m = FN_RE.match(os.path.basename(f))
    if not m: return None
    y, mo, d = map(int, m.groups()[:3])
    return date(y, mo, d)

def find_files_in_range(path: str, d0: date, d1: date) -> List[str]:
    pats = [os.path.join(path, "iceh_inst.*.nc"),
            os.path.join(path, "*", "iceh_inst.*.nc")]
    files = []
    for p in pats:
        files.extend(glob.glob(p))
    keep = []
    for f in files:
        fd = file_date(f)
        if fd and (d0 <= fd <= d1):
            keep.append(f)
    keep.sort()
    return keep

def read_first(ds: xr.Dataset, names: List[str]) -> Optional[xr.DataArray]:
    """Return first existing var as float DataArray (time squeezed, NaN for non-finite)."""
    for n in names:
        if n in ds:
            da = ds[n]
            # squeeze time if present and length==1
            if "time" in da.dims and da.sizes.get("time", None) == 1:
                da = da.isel(time=0, drop=True)
            # ensure float and mask non-finite
            da = da.astype("float64")
            da = da.where(np.isfinite(da))
            return da
    return None

def nan_max(arr: np.ndarray) -> float:
    """NaN-safe absolute max; returns np.nan if all-NaN or empty."""
    a = np.asarray(arr)
    a = a[np.isfinite(a)]
    return float(np.max(np.abs(a))) if a.size else float("nan")

def nan_p99(arr: np.ndarray) -> float:
    """NaN-safe 99th percentile of absolute; returns np.nan if all-NaN or empty."""
    a = np.asarray(arr)
    a = a[np.isfinite(a)]
    return float(np.percentile(np.abs(a), 99)) if a.size else float("nan")

def analyze_file_box(nc: str) -> Optional[Dict[str, float]]:
    try:
        ds = xr.open_dataset(nc, decode_cf=True, mask_and_scale=True)
    except Exception as e:
        print(f"[WARN] Failed to open {os.path.basename(nc)}: {e}")
        return None

    # Required face vars (E/N normals & tangentials) + masks (optional/fallbacks handled)
    uE = read_first(ds, ["uvelE_1","UVEL_E","uE"])
    vE = read_first(ds, ["vvelE_1","VVEL_E","vE"])
    uN = read_first(ds, ["uvelN_1","UVEL_N","uN"])
    vN = read_first(ds, ["vvelN_1","VVEL_N","vN"])

    # Masks may be missing or NaN; we'll tolerate that
    emask = read_first(ds, ["emask","EMASK"])
    nmask = read_first(ds, ["nmask","NMASK"])

    # All required velocity fields must exist
    if any(x is None for x in (uE, vE, uN, vN)):
        print(f"[WARN] Missing required face velocities in {os.path.basename(nc)}")
        ds.close()
        return None

    # Ensure shapes are (nj, ni)
    try:
        nj, ni = uE.shape
    except Exception:
        # attempt to align by dropping any extra dims except (nj, ni)
        for name, da in [("uE", uE), ("vE", vE), ("uN", uN), ("vN", vN)]:
            while da.ndim > 2:
                da = da.isel({da.dims[0]: 0}, drop=True)
            if name == "uE": uE = da
            if name == "vE": vE = da
            if name == "uN": uN = da
            if name == "vN": vN = da
        nj, ni = uE.shape

    jj, ii = np.indices((nj, ni))

    # Ocean masks for faces (robust to NaN): True where mask >= 0.5; if mask missing/unusable, fallback to "finite velocity"
    if emask is not None and np.isfinite(emask.values).any():
        oceanE = (emask >= 0.5).fillna(False).values
    else:
        oceanE = np.isfinite(vE.values) | np.isfinite(uE.values)

    if nmask is not None and np.isfinite(nmask.values).any():
        oceanN = (nmask >= 0.5).fillna(False).values
    else:
        oceanN = np.isfinite(vN.values) | np.isfinite(uN.values)

    # Box-specific “coast”: outermost faces only, intersected with ocean mask.
    # If intersect kills everything (e.g., masks are all-NaN), fall back to edges alone.
    coastE = (((ii == 0) | (ii == ni-1)) & oceanE)
    coastN = (((jj == 0) | (jj == nj-1)) & oceanN)
    if coastE.sum() == 0:
        coastE = ((ii == 0) | (ii == ni-1))
    if coastN.sum() == 0:
        coastN = ((jj == 0) | (jj == nj-1))

    # Interior (ocean) faces away from edges for reference scales
    E_int = (oceanE & (ii > 0) & (ii < ni-1))
    N_int = (oceanN & (jj > 0) & (jj < nj-1))

    # Reference scales (like-for-like, NaN-safe)
    uE_ref = nan_p99(uE.values[E_int])   # E normals
    vN_ref = nan_p99(vN.values[N_int])   # N normals
    vE_ref = nan_p99(vE.values[E_int])   # E tangential
    uN_ref = nan_p99(uN.values[N_int])   # N tangential

    # Coast selections (NaN elsewhere)
    uE_coast = np.where(coastE, uE.values, np.nan)   # normal at E
    vN_coast = np.where(coastN, vN.values, np.nan)   # normal at N
    vE_coast = np.where(coastE, vE.values, np.nan)   # tangential at E
    uN_coast = np.where(coastN, uN.values, np.nan)   # tangential at N

    # Free-slip “normal derivative” (index units) using nearest interior neighbor
    dv_dn = np.full_like(vE.values, np.nan, dtype=float)
    du_dn = np.full_like(uN.values, np.nan, dtype=float)

    # E: west edge (i=0) → compare to i=1; east edge (i=ni-1) → compare to i=ni-2
    west = coastE & (ii == 0)
    east = coastE & (ii == ni-1)
    if ni >= 2:
        dv_dn[west] = vE.values[west] - vE.values[jj[west], 1]
        dv_dn[east] = vE.values[east] - vE.values[jj[east], ni-2]

    # N: south edge (j=0) → compare to j=1; north edge (j=nj-1) → compare to j=nj-2
    south = coastN & (jj == 0)
    north = coastN & (jj == nj-1)
    if nj >= 2:
        du_dn[south] = uN.values[south] - uN.values[1, ii[south]]
        du_dn[north] = uN.values[north] - uN.values[nj-2, ii[north]]

    out = dict(
        cntE=int(np.count_nonzero(coastE)),
        cntN=int(np.count_nonzero(coastN)),
        uE_ref=uE_ref, vN_ref=vN_ref, vE_ref=vE_ref, uN_ref=uN_ref,
        max_un_E=nan_max(uE_coast),
        max_vn_N=nan_max(vN_coast),
        max_vt_E=nan_max(vE_coast),
        max_ut_N=nan_max(uN_coast),
        max_dv_dn=nan_max(dv_dn),
        max_du_dn=nan_max(du_dn),
    )
    ds.close()
    return out

def main():
    ap = argparse.ArgumentParser(description="NaN-aware box BC checker (all 4 edges coastline).")
    ap.add_argument("--start", required=True, help="YYYY-MM-DD")
    ap.add_argument("--end",   required=True, help="YYYY-MM-DD")
    ap.add_argument("--exp", action="append", default=[],
                    help="LABEL:PATH (history dir). Repeatable. LABEL in {free,noslip,noratio,any}.")
    args = ap.parse_args()

    d0 = parse_datestr(args.start)
    d1 = parse_datestr(args.end)
    if d1 < d0:
        print("ERROR: end date before start date"); sys.exit(2)

    if not args.exp:
        print("Provide at least one --exp LABEL:PATH"); sys.exit(2)

    for item in args.exp:
        if ":" not in item:
            print(f"Bad --exp '{item}' (need LABEL:PATH)"); continue
        label, path = item.split(":",1)
        path = path.strip(); label = label.strip().upper()
        files = find_files_in_range(path, d0, d1)
        usable = []
        per = []
        for f in files:
            res = analyze_file_box(f)
            if res is not None:
                usable.append(f); per.append(res)

        if not files:
            print(f"=== {label} () ==="); print(f"No files in range under {path}\n"); continue
        if not per:
            print(f"=== {label} () ==="); print(f"Files in range: {len(files)} | usable: 0\n"); continue

        # average of per-file maxima (NaN-safe)
        def avg(key):
            vals = np.array([r[key] for r in per], float)
            vals = vals[np.isfinite(vals)]
            return float(vals.mean()) if vals.size else float("nan")

        cntE_min = min(r["cntE"] for r in per); cntE_max = max(r["cntE"] for r in per)
        cntN_min = min(r["cntN"] for r in per); cntN_max = max(r["cntN"] for r in per)

        print(f"=== {label} () ===")
        print(f"Files in range: {len(files)} | usable: {len(usable)} | coastal faces E: {cntE_min}-{cntE_max}, N: {cntN_min}-{cntN_max}")
        print("Avg interior scales (99th pct):")
        print(f"  normals:  uE(E)={avg('uE_ref'):0.3e}, vN(N)={avg('vN_ref'):0.3e}")
        print(f"  tang/der: vE(E)={avg('vE_ref'):0.3e}, uN(N)={avg('uN_ref'):0.3e}")

        print("Normals ~ 0 (both BCs): (average of per-file maxima)")
        print(f"  avg max |u_n| on E coast   : {avg('max_un_E'):0.3e}")
        print(f"  avg max |v_n| on N coast   : {avg('max_vn_N'):0.3e}")

        print("Tangentials on coast (no-slip ~0; free-slip unconstrained):")
        vtE = avg('max_vt_E'); utN = avg('max_ut_N')
        vE_ref = avg('vE_ref'); uN_ref = avg('uN_ref')
        rel_vtE = (vtE/vE_ref) if np.isfinite(vtE) and np.isfinite(vE_ref) and vE_ref>0 else float('nan')
        rel_utN = (utN/uN_ref) if np.isfinite(utN) and np.isfinite(uN_ref) and uN_ref>0 else float('nan')
        print(f"  avg max |v_tan| on E coast : {vtE:0.3e}  (rel={rel_vtE:0.3e})")
        print(f"  avg max |u_tan| on N coast : {utN:0.3e}  (rel={rel_utN:0.3e})")

        print("Normal derivative of tangential (free-slip ~0; index units):")
        print(f"  avg max |∂v/∂n| on E coast : {avg('max_dv_dn'):0.3e}")
        print(f"  avg max |∂u/∂n| on N coast : {avg('max_du_dn'):0.3e}")

        # Soft pass/fail (tuned for box) — NaN-aware and skip if no coastal faces
        tol_rel = 1e-3
        abs_floor_deriv = 1e-7

        def rel_ok(val, ref):
            if not np.isfinite(val): return False
            if not np.isfinite(ref) or ref <= 0: return True
            return val <= tol_rel*ref

        def deriv_ok(val, ref):
            if not np.isfinite(val): return False
            if val <= abs_floor_deriv: return True
            if not np.isfinite(ref) or ref <= 0: return True
            return val <= 1e-2*ref

        normals_ok = rel_ok(avg('max_un_E'), avg('uE_ref')) and rel_ok(avg('max_vn_N'), avg('vN_ref'))
        noslip_ok  = rel_ok(vtE, vE_ref) and rel_ok(utN, uN_ref)
        freeslip_ok = deriv_ok(avg('max_dv_dn'), vE_ref) and deriv_ok(avg('max_du_dn'), uN_ref)

        # If no coastal faces detected at all, don't judge
        if cntE_max == 0 and cntN_max == 0:
            print("Verdict: no coastal faces detected (skipping PASS/FAIL)\n")
            continue

        lab = label.lower()
        if "noslip" in lab or "no-slip" in lab:
            print(f"Verdict (no-slip): normals={'OK' if normals_ok else 'FAIL'}, tangentials={'OK' if noslip_ok else 'FAIL'}\n")
        elif "free" in lab or "noratio" in lab:
            print(f"Verdict (free-slip): normals={'OK' if normals_ok else 'FAIL'}, d(tangent)/dn={'OK' if freeslip_ok else 'FAIL'}  (tangentials may be O(1) of interior)\n")
        else:
            print(f"Verdict: normals={'OK' if normals_ok else 'FAIL'}; noslip tang={'OK' if noslip_ok else 'FAIL'}; freeslip d/dn={'OK' if freeslip_ok else 'FAIL'}\n")

if __name__ == "__main__":
    main()
