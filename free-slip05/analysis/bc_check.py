#!/usr/bin/env python3
"""
Boundary-condition checks for CICE instantaneous files (C-grid), over a time range.

For each experiment directory, this script:
  • finds all iceh_inst.YYYY-MM-DD-SSSSS.nc between --start and --end (inclusive),
  • for each file: computes
      - coastal masks from T mask,
      - normals on coast (u on E faces, v on N faces),
      - tangentials on coast (v on E, u on N),
      - free-slip “normal derivative of tangential” in index units:
          E: v(j,i) - v(j, i±1 toward interior)
          N: u(j,i) - u(j±1, i toward interior)
      - interior reference scales (99th percentile away from coast) for v on E, u on N,
    with full NaN handling,
  • aggregates over time by averaging each per-file maximum (average of maxima).

Usage:
  python bc_check.py \
      --start 2005-01-01 --end 2005-01-03 \
      --exp free:/path/to/BC_free/history \
      --exp noratio:/path/to/BC_free_noratio/history \
      --exp noslip:/path/to/BC_noslip/history

Notes
  • Variables expected: tmask, emask, nmask, uvelE_1, vvelE_1, uvelN_1, vvelN_1.
    (Name fallbacks supported: UVEL_E/VVEL_E/uE/vE, UVEL_N/VVEL_N/uN/vN).
  • Derivatives are in *index* units (grid spacing not used unless you adapt dx,dy).
"""

import argparse, glob, os, re
from datetime import datetime, timedelta
from typing import Optional, List, Tuple, Dict
import numpy as np
from netCDF4 import Dataset

FN_RE = re.compile(r"iceh_inst\.(\d{4})-(\d{2})-(\d{2})-(\d{5})\.nc$")

def parse_file_dt(path: str) -> Optional[datetime]:
    m = FN_RE.search(os.path.basename(path))
    if not m: return None
    y, mo, d, sssss = map(int, m.groups())
    return datetime(y, mo, d) + timedelta(seconds=sssss)

def parse_user_time(s: str) -> datetime:
    s = s.strip()
    for fmt in ("%Y-%m-%d", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d-%H%M%S"):
        try: return datetime.strptime(s, fmt)
        except ValueError: pass
    m = re.match(r"^(\d{4}-\d{2}-\d{2})-(\d{5})$", s)
    if m:
        d = datetime.strptime(m.group(1), "%Y-%m-%d")
        return d + timedelta(seconds=int(m.group(2)))
    raise ValueError(f"Unrecognized time format: {s}")

def read_var(ds: Dataset, names: List[str]) -> Optional[np.ndarray]:
    for name in names:
        if name in ds.variables:
            v = ds.variables[name]
            a = v[:]
            a = a.filled(np.nan) if np.ma.isMaskedArray(a) else np.asarray(a, float)
            for attr in ("_FillValue","missing_value"):
                if attr in v.ncattrs():
                    try:
                        fv = float(getattr(v, attr))
                        a = np.where(np.isfinite(a) & (a == fv), np.nan, a)
                    except Exception:
                        pass
            return a.astype(float, copy=False)
    return None

def squeeze_time(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a)
    return a[0] if (a.ndim == 3 and a.shape[0] == 1) else a

def xor(a,b): return np.logical_xor(a,b)

def safe_percentile(a: np.ndarray, p: float, default=np.nan) -> float:
    a = np.asarray(a); a = a[np.isfinite(a)]
    return default if a.size==0 else float(np.percentile(np.abs(a), p))

def safe_max(a: np.ndarray, default=np.nan) -> float:
    a = np.asarray(a); a = a[np.isfinite(a)]
    return default if a.size==0 else float(np.max(np.abs(a)))

# ---------- NEW: robust coast detection (internal + perimeter) ----------
def build_coast_masks(tmask: Optional[np.ndarray],
                      emask: np.ndarray,
                      nmask: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Returns coastE, coastN (boolean) with the same (nj,ni) shape as *_mask vars.
    A face is 'coast' if:
      - it’s an ocean face AND sits on the domain perimeter, OR
      - it’s an ocean face and its adjacent T cells across the face differ (if tmask given).
    """
    Eo = emask > 0.5
    No = nmask > 0.5
    nj, ni = Eo.shape
    coastE = np.zeros((nj, ni), bool)
    coastN = np.zeros((nj, ni), bool)

    # Perimeter-as-coast for an ocean face
    coastE[:, 0]  |= Eo[:, 0]
    coastE[:, -1] |= Eo[:, -1]
    coastN[0, :]  |= No[0, :]
    coastN[-1, :] |= No[-1, :]

    # Internal coasts from T mask if available
    if tmask is not None:
        T = tmask > 0.0
        coastE[:, :-1] |= Eo[:, :-1] & xor(T[:, :-1], T[:, 1:])
        coastN[:-1, :] |= No[:-1, :] & xor(T[:-1, :], T[1:, :])

    # Also mark faces adjacent to a land face in the normal direction (emask/nmask)
    # E: if left or right neighbor E-face is land (0/NaN)
    left_is_land  = np.zeros_like(Eo, bool); left_is_land[:, 1:]  = ~ (emask[:, :-1] > 0.5)
    right_is_land = np.zeros_like(Eo, bool); right_is_land[:, :-1] = ~ (emask[:, 1:]  > 0.5)
    coastE |= Eo & (left_is_land | right_is_land)

    # N: if below or above neighbor N-face is land
    below_is_land = np.zeros_like(No, bool); below_is_land[1:, :]  = ~ (nmask[:-1, :] > 0.5)
    above_is_land = np.zeros_like(No, bool); above_is_land[:-1, :] = ~ (nmask[1:,  :] > 0.5)
    coastN |= No & (below_is_land | above_is_land)

    return coastE, coastN

# ---------- NEW: pick interior neighbor for free-slip derivative ----------
def interior_offset_E(tmask: Optional[np.ndarray], emask: np.ndarray) -> np.ndarray:
    nj, ni = emask.shape
    off = np.zeros((nj, ni), np.int8)
    # Prefer direction where adjacent E-face is ocean; fallback to edges; use tmask if present.
    if tmask is not None:
        T = tmask > 0.0
        left_ocean  = np.zeros_like(T, bool);  left_ocean[:, 1:]  =  T[:, 1:] & (~T[:, :-1])
        right_ocean = np.zeros_like(T, bool);  right_ocean[:, :-1] = (~T[:, :-1]) & T[:, 1:]
        off[left_ocean]  = -1   # ocean on left -> interior to left
        off[right_ocean] = +1   # ocean on right -> interior to right
    # Use emask around the face
    left_is_ocean  = np.zeros_like(emask, bool); left_is_ocean[:, 1:]  = emask[:, :-1] > 0.5
    right_is_ocean = np.zeros_like(emask, bool); right_is_ocean[:, :-1] = emask[:, 1:]  > 0.5
    off[(right_is_ocean) & (~left_is_ocean)] = +1
    off[(left_is_ocean)  & (~right_is_ocean)] = -1
    # Perimeter fallback
    off[:, 0]  = np.where(emask[:, 0]  > 0.5, +1, 0)
    off[:, -1] = np.where(emask[:, -1] > 0.5, -1, 0)
    return off

def interior_offset_N(tmask: Optional[np.ndarray], nmask: np.ndarray) -> np.ndarray:
    nj, ni = nmask.shape
    off = np.zeros((nj, ni), np.int8)
    if tmask is not None:
        T = tmask > 0.0
        below_ocean = np.zeros_like(T, bool); below_ocean[1:, :]  =  T[1:, :] & (~T[:-1, :])
        above_ocean = np.zeros_like(T, bool); above_ocean[:-1, :] = (~T[:-1, :]) & T[1:, :]
        off[below_ocean] = -1   # ocean below -> interior downward
        off[above_ocean] = +1   # ocean above -> interior upward
    # Use nmask around the face
    below_is_ocean = np.zeros_like(nmask, bool); below_is_ocean[1:, :]  = nmask[:-1, :] > 0.5
    above_is_ocean = np.zeros_like(nmask, bool); above_is_ocean[:-1, :] = nmask[1:,  :] > 0.5
    off[(above_is_ocean) & (~below_is_ocean)] = +1
    off[(below_is_ocean) & (~above_is_ocean)] = -1
    # Perimeter fallback
    off[0, :]  = np.where(nmask[0, :]  > 0.5, +1, 0)
    off[-1, :] = np.where(nmask[-1, :] > 0.5, -1, 0)
    return off

def list_inst_files_between(hist_dir: str, start_dt: datetime, end_dt: datetime) -> List[str]:
    pats = [os.path.join(hist_dir, "iceh_inst.*.nc"),
            os.path.join(hist_dir, "*", "iceh_inst.*.nc")]
    files = []
    for p in pats: files.extend(glob.glob(p))
    keep = []
    for f in files:
        dt = parse_file_dt(f)
        if dt and (start_dt <= dt <= end_dt):
            keep.append((dt, f))
    keep.sort(key=lambda x: x[0])
    return [f for _, f in keep]

def analyze_one_file(nc_path: str) -> dict:
    ds = Dataset(nc_path, "r")
    tmask = read_var(ds, ["tmask","TMSK","tmask_1"])
    emask = read_var(ds, ["emask","EMASK"])
    nmask = read_var(ds, ["nmask","NMASK"])
    uE = read_var(ds, ["uvelE_1","UVEL_E","uE"])
    vE = read_var(ds, ["vvelE_1","VVEL_E","vE"])
    uN = read_var(ds, ["uvelN_1","UVEL_N","uN"])
    vN = read_var(ds, ["vvelN_1","VVEL_N","vN"])

    need = (emask, nmask, uE, vE, uN, vN)
    if any(x is None for x in need):
        ds.close(); return {"ok": False, "reason": "missing variables"}

    # time squeeze
    if tmask is not None: tmask = squeeze_time(tmask).astype(float)
    emask = squeeze_time(emask).astype(float)
    nmask = squeeze_time(nmask).astype(float)
    uE = squeeze_time(uE).astype(float)
    vE = squeeze_time(vE).astype(float)
    uN = squeeze_time(uN).astype(float)
    vN = squeeze_time(vN).astype(float)

    coastE, coastN = build_coast_masks(tmask, emask, nmask)

    # interior (away from coasts) for reference scales
    E_notcoast = (emask > 0.5) & (~coastE)
    N_notcoast = (nmask > 0.5) & (~coastN)
    vE_ref = safe_percentile(vE[E_notcoast], 99)
    uN_ref = safe_percentile(uN[N_notcoast], 99)

    # normals/tangentials restricted to coast
    u_n_E = np.where(coastE, uE, np.nan)
    v_n_N = np.where(coastN, vN, np.nan)
    v_tan_E = np.where(coastE, vE, np.nan)
    u_tan_N = np.where(coastN, uN, np.nan)

    # free-slip derivative in index units
    offE = interior_offset_E(tmask, emask)
    offN = interior_offset_N(tmask, nmask)
    dv_dn = np.full_like(vE, np.nan, float)
    du_dn = np.full_like(uN, np.nan, float)
    jj, ii = np.indices(vE.shape)
    ii_nei = np.clip(ii + offE, 0, vE.shape[1]-1)
    jj_nei = np.clip(jj + offN, 0, uN.shape[0]-1)
    validE = coastE & (offE != 0)
    validN = coastN & (offN != 0)
    dv_dn[validE] = vE[validE] - vE[jj[validE], ii_nei[validE]]
    du_dn[validN] = uN[validN] - uN[jj_nei[validN], ii[validN]]

    out = {
        "ok": True,
        "cnt_E": int(np.count_nonzero(coastE)),
        "cnt_N": int(np.count_nonzero(coastN)),
        "vE_ref": vE_ref, "uN_ref": uN_ref,
        "max_un_E": safe_max(u_n_E), "max_vn_N": safe_max(v_n_N),
        "max_vt_E": safe_max(v_tan_E), "max_ut_N": safe_max(u_tan_N),
        "max_dv_dn": safe_max(dv_dn), "max_du_dn": safe_max(du_dn),
    }
    ds.close()
    return out

def nanmean(xs: List[float]) -> float:
    a = np.array(xs, float); a = a[np.isfinite(a)]
    return np.nan if a.size==0 else float(np.mean(a))

def summarize_over_time(per_file: List[dict]):
    ok = [d for d in per_file if d.get("ok", False)]
    if not ok: return None
    agg = {}
    for k in ("max_un_E","max_vn_N","max_vt_E","max_ut_N","max_dv_dn","max_du_dn","vE_ref","uN_ref"):
        agg[k] = nanmean([d[k] for d in ok])
    agg["cnt_E_min"] = min(d["cnt_E"] for d in ok)
    agg["cnt_E_max"] = max(d["cnt_E"] for d in ok)
    agg["cnt_N_min"] = min(d["cnt_N"] for d in ok)
    agg["cnt_N_max"] = max(d["cnt_N"] for d in ok)
    agg["nfiles_ok"] = len(ok); agg["nfiles_total"] = len(per_file)
    return agg

def main():
    ap = argparse.ArgumentParser(description="Check CICE BCs over a time span (impermeability, no-/free-slip).")
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--exp", action="append", default=[],
                    help="LABEL:PATH (LABEL in {free,noratio,noslip,any})")
    args = ap.parse_args()

    start_dt = parse_user_time(args.start)
    end_dt   = parse_user_time(args.end)
    if end_dt < start_dt:
        print("End time precedes start time."); return 2
    if not args.exp:
        print("Provide at least one --exp LABEL:PATH"); return 2

    label_to_expect = {
        "noslip":"noslip","no-slip":"noslip",
        "free":"free","free-slip":"free",
        "noratio":"noratio","free-slip-noratio":"noratio",
    }

    exps: Dict[str, Tuple[str,str]] = {}
    for item in args.exp:
        if ":" not in item:
            print(f"Bad --exp '{item}'. Expected LABEL:PATH"); return 2
        label, path = item.split(":",1)
        label = label.strip(); path = path.strip()
        if not os.path.isdir(path):
            print(f"[{label}] Not a directory: {path}"); continue
        exps[label] = (path, label_to_expect.get(label.lower(),"any"))
    if not exps:
        print("No valid experiments found."); return 1

    for label, (path, expect) in exps.items():
        files = list_inst_files_between(path, start_dt, end_dt)
        print(f"=== {label.upper()} ({os.path.basename(path)}) ===")
        if not files:
            print(f"No files in range under: {path}\n"); continue
        per_file = [analyze_one_file(f) for f in files]
        agg = summarize_over_time(per_file)
        if not agg:
            print("No usable files (missing variables?).\n"); continue

        print(f"Files in range: {len(files)} | usable: {agg['nfiles_ok']} "
              f"| coastal faces E: {agg['cnt_E_min']}-{agg['cnt_E_max']}, "
              f"N: {agg['cnt_N_min']}-{agg['cnt_N_max']}")
        vE_ref, uN_ref = agg["vE_ref"], agg["uN_ref"]
        rel_vt_E = (agg["max_vt_E"]/vE_ref) if (np.isfinite(vE_ref) and vE_ref>0) else np.nan
        rel_ut_N = (agg["max_ut_N"]/uN_ref) if (np.isfinite(uN_ref) and uN_ref>0) else np.nan

        print(f"Avg interior scales (99th pct): V_E={vE_ref:0.3e} (v on E), U_N={uN_ref:0.3e} (u on N)")
        print("Normals ~ 0 (both BCs): (average of per-file maxima)")
        print(f"  avg max |u_n| on E coast   : {agg['max_un_E']:0.3e}")
        print(f"  avg max |v_n| on N coast   : {agg['max_vn_N']:0.3e}")
        print("Tangentials on coast (no-slip ~0; free-slip unconstrained):")
        print(f"  avg max |v_tan| on E coast : {agg['max_vt_E']:0.3e}  (rel={rel_vt_E:0.3e})")
        print(f"  avg max |u_tan| on N coast : {agg['max_ut_N']:0.3e}  (rel={rel_ut_N:0.3e})")
        print("Normal derivative of tangential (free-slip ~0; index units):")
        print(f"  avg max |∂v/∂n| on E coast : {agg['max_dv_dn']:0.3e}")
        print(f"  avg max |∂u/∂n| on N coast : {agg['max_du_dn']:0.3e}")

        tol_rel = 1e-3
        def rel_ok(val, scale):
            if not np.isfinite(val): return False
            if not np.isfinite(scale) or scale <= 0: return True
            return val <= tol_rel*scale

        normals_ok  = rel_ok(agg["max_un_E"], max(vE_ref, uN_ref)) and rel_ok(agg["max_vn_N"], max(vE_ref, uN_ref))
        noslip_ok   = rel_ok(agg["max_vt_E"], vE_ref) and rel_ok(agg["max_ut_N"], uN_ref)
        freeslip_ok = rel_ok(agg["max_dv_dn"], vE_ref) and rel_ok(agg["max_du_dn"], uN_ref)

        if expect == "noslip":
            print(f"Verdict (no-slip): normals={'OK' if normals_ok else 'FAIL'}, "
                  f"tangentials={'OK' if noslip_ok else 'FAIL'}\n")
        elif expect in ("free","noratio"):
            print(f"Verdict (free-slip): normals={'OK' if normals_ok else 'FAIL'}, "
                  f"d(tangent)/dn={'OK' if freeslip_ok else 'FAIL'}  "
                  f"(tangentials may be O(1) of interior)\n")
        else:
            print(f"Verdict: normals={'OK' if normals_ok else 'FAIL'}; "
                  f"no-slip tangentials={'OK' if noslip_ok else 'FAIL'}; "
                  f"free-slip d(tan)/dn={'OK' if freeslip_ok else 'FAIL'}\n")

if __name__ == "__main__":
    raise SystemExit(main())


