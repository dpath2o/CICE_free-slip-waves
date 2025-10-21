#!/usr/bin/env python3

"""
Boundary-condition checks for CICE instantaneous files (C-grid).
Verifies impermeability (normals), no-slip (tangentials), and free-slip (normal derivative of tangential).
Works with variables commonly present in iceh_inst.*.nc outputs (see README in code).

Usage examples:
  python bc_check.py --exp free:/path/to/BC_free/history \
                     --exp noratio:/path/to/BC_free_noratio/history \
                     --exp noslip:/path/to/BC_noslip/history

By default, the newest iceh_inst.*.nc in each experiment directory is selected.

Outputs compact text summaries, with NaN-safe handling.
"""

import argparse
import glob
import os
import sys
from typing import Dict, Tuple, Optional

import numpy as np
from netCDF4 import Dataset

# -------------------------------
# Helpers
# -------------------------------

def find_latest_inst_file(hist_dir: str) -> Optional[str]:
    pats = [os.path.join(hist_dir, "iceh_inst.*.nc"),
            os.path.join(hist_dir, "*", "iceh_inst.*.nc")]
    files = []
    for pat in pats:
        files.extend(glob.glob(pat))
    if not files:
        return None
    files.sort()
    return files[-1]

def read_var(ds: Dataset, names):
    """Try a list of candidate variable names and return the first that exists; else None."""
    for name in names:
        if name in ds.variables:
            return ds.variables[name][:]
    return None

def valid_mask(arr):
    """Finite mask for array (True where finite)."""
    return np.isfinite(arr)

def xor(a, b):
    return np.logical_xor(a, b)

def safe_percentile(a: np.ndarray, p: float, default: float = 0.0) -> float:
    a = np.asarray(a)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return default
    return float(np.percentile(np.abs(a), p))

def safe_max(a: np.ndarray, default: float = np.nan) -> float:
    a = np.asarray(a)
    a = a[np.isfinite(a)]
    if a.size == 0:
        return default
    return float(np.max(np.abs(a)))

def build_coast_masks_from_tmask(tmask: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build coastal face masks from T-grid mask.
    E face (j,i) is the vertical interface between T(i,j) and T(i+1,j).
    N face (j,i) is the horizontal interface between T(i,j) and T(i,j+1).

    Returns:
      coastE (nj,ni) with last col False,
      coastN (nj,ni) with last row False.
    """
    ocean = np.asarray(tmask) > 0.0  # boolean
    nj, ni = ocean.shape

    coastE = np.zeros((nj, ni), dtype=bool)
    coastN = np.zeros((nj, ni), dtype=bool)

    # E faces exist for i = 0..ni-2
    coastE[:, :-1] = xor(ocean[:, :-1], ocean[:, 1:])

    # N faces exist for j = 0..nj-2
    coastN[:-1, :] = xor(ocean[:-1, :], ocean[1:, :])

    return coastE, coastN

def neighbor_interior_index_E(tmask: np.ndarray) -> np.ndarray:
    """
    For each coastal E face (j,i), choose the adjacent interior E face index offset in x:
      if ocean is on left (tmask[j,i]==1, tmask[j,i+1]==0) -> offset = -1 (i-1)
      if ocean is on right (tmask[j,i]==0, tmask[j,i+1]==1) -> offset = +1 (i+1)
    Non-coastal or edges with no interior neighbor get offset 0.

    Returns an int array shape (nj,ni) with values in {-1,0,+1}.
    """
    ocean = tmask > 0.0
    nj, ni = ocean.shape
    off = np.zeros((nj, ni), dtype=np.int8)

    left_ocean  = np.zeros_like(ocean, dtype=bool)
    right_ocean = np.zeros_like(ocean, dtype=bool)
    left_ocean[:, :-1]  = ocean[:, :-1] & (~ocean[:, 1:])
    right_ocean[:, :-1] = (~ocean[:, :-1]) & ocean[:, 1:]

    off[left_ocean]  = -1
    off[right_ocean] = +1

    # Remove places where neighbor index would be out of bounds
    off[:, 0] = 0
    off[:, -1] = 0
    return off

def neighbor_interior_index_N(tmask: np.ndarray) -> np.ndarray:
    """
    For each coastal N face (j,i), choose the adjacent interior N face index offset in y:
      if ocean is below (tmask[j,i]==1, tmask[j+1,i]==0) -> offset = -1 (j-1)
      if ocean is above (tmask[j,i]==0, tmask[j+1,i]==1) -> offset = +1 (j+1)
    Non-coastal or edges with no interior neighbor get offset 0.

    Returns an int array shape (nj,ni) with values in {-1,0,+1}.
    """
    ocean = tmask > 0.0
    nj, ni = ocean.shape
    off = np.zeros((nj, ni), dtype=np.int8)

    below_ocean = np.zeros_like(ocean, dtype=bool)
    above_ocean = np.zeros_like(ocean, dtype=bool)
    below_ocean[:-1, :] = ocean[:-1, :] & (~ocean[1:, :])
    above_ocean[:-1, :] = (~ocean[:-1, :]) & ocean[1:, :]

    off[below_ocean] = -1
    off[above_ocean] = +1

    off[0, :]  = 0
    off[-1, :] = 0
    return off

def analyze_one_file(nc_path: str, label: str, expect: str):
    """
    expect: one of {"free", "noratio", "noslip"} (used only for PASS/FAIL hints)
    """
    ds = Dataset(nc_path, "r")

    # Pull masks and velocities (support a few name variants)
    tmask = read_var(ds, ["tmask", "TMSK", "tmask_1"])
    emask = read_var(ds, ["emask", "EMASK"])
    nmask = read_var(ds, ["nmask", "NMASK"])

    uE = read_var(ds, ["uvelE_1", "UVEL_E", "uE"])
    vE = read_var(ds, ["vvelE_1", "VVEL_E", "vE"])
    uN = read_var(ds, ["uvelN_1", "UVEL_N", "uN"])
    vN = read_var(ds, ["vvelN_1", "VVEL_N", "vN"])

    # Optional metrics (not required; derivatives are in index units)
    # dxu = read_var(ds, ["dxu", "DXU"]); dyu = read_var(ds, ["dyu", "DYU"])

    if any(x is None for x in (tmask, emask, nmask, uE, vE, uN, vN)):
        print(f"[{label}] ERROR: required variables missing in {os.path.basename(nc_path)}")
        ds.close()
        return

    # Squeeze time if present
    def squeeze2(a):
        a = np.asarray(a)
        if a.ndim == 3:  # time,nj,ni
            a = a[0]
        return a

    tmask = squeeze2(tmask).astype(float)
    emask = squeeze2(emask).astype(float)
    nmask = squeeze2(nmask).astype(float)
    uE = squeeze2(uE).astype(float)
    vE = squeeze2(vE).astype(float)
    uN = squeeze2(uN).astype(float)
    vN = squeeze2(vN).astype(float)

    nj, ni = tmask.shape

    # Build coastal masks from T mask
    coastE, coastN = build_coast_masks_from_tmask(tmask)

    # Ocean masks on faces
    oceanE = emask > 0.5
    oceanN = nmask > 0.5

    # Coastal *ocean* faces
    cE = coastE & oceanE
    cN = coastN & oceanN

    # Interior (ocean) faces one cell away from coast (for reference scale)
    # E interior: ocean & not coastal & both neighbors ocean to avoid edge effects
    E_notcoast = oceanE & (~coastE)
    N_notcoast = oceanN & (~coastN)

    # Interior reference scales (99th percentile) for tangentials
    vE_ref = safe_percentile(vE[E_notcoast], 99, default=0.0)
    uN_ref = safe_percentile(uN[N_notcoast], 99, default=0.0)

    # -------- Tests --------
    # 1) Impermeability: normal components at coast ~ 0
    u_n_E = np.where(cE, uE, np.nan)
    v_n_N = np.where(cN, vN, np.nan)

    # 2) Tangentials at coast
    v_tan_E = np.where(cE, vE, np.nan)
    u_tan_N = np.where(cN, uN, np.nan)

    # 3) Free-slip derivative: d(v)/dn on E, d(u)/dn on N (index units)
    offE = neighbor_interior_index_E(tmask)  # {-1,0,+1}
    offN = neighbor_interior_index_N(tmask)  # {-1,0,+1}

    dv_dn = np.full_like(vE, np.nan, dtype=float)
    du_dn = np.full_like(uN, np.nan, dtype=float)

    # Compute using nearest interior neighbor in x for E
    jj, ii = np.indices(vE.shape)
    ii_nei = ii + offE
    validE = cE & (offE != 0)
    ii_nei = np.clip(ii_nei, 0, ni - 1)
    dv_dn[validE] = vE[validE] - vE[jj[validE], ii_nei[validE]]

    # Compute using nearest interior neighbor in y for N
    jj_nei = jj + offN
    validN = cN & (offN != 0)
    jj_nei = np.clip(jj_nei, 0, nj - 1)
    du_dn[validN] = uN[validN] - uN[jj_nei[validN], ii[validN]]

    # Summaries (NaN-safe)
    cnt_E = int(np.count_nonzero(cE))
    cnt_N = int(np.count_nonzero(cN))

    max_un_E = safe_max(u_n_E, default=np.nan)
    max_vn_N = safe_max(v_n_N, default=np.nan)
    max_vt_E = safe_max(v_tan_E, default=np.nan)
    max_ut_N = safe_max(u_tan_N, default=np.nan)
    max_dv_dn = safe_max(dv_dn, default=np.nan)
    max_du_dn = safe_max(du_dn, default=np.nan)

    # Relative magnitudes vs interior scales (avoid div by 0)
    rel_vt_E = (max_vt_E / vE_ref) if vE_ref > 0 else np.nan
    rel_ut_N = (max_ut_N / uN_ref) if uN_ref > 0 else np.nan

    # Print
    print(f"=== {label.upper()} ({os.path.basename(nc_path)}) ===")
    print(f"Coast E faces count: {cnt_E}; Coast N faces count: {cnt_N}")
    print(f"Interior scales (99th pct): V_E={vE_ref:0.3e} (v on E), U_N={uN_ref:0.3e} (u on N)")
    print("Normals ~ 0 (both BCs):")
    print(f"  max |u_n| on E coast   : {max_un_E:0.3e}")
    print(f"  max |v_n| on N coast   : {max_vn_N:0.3e}")
    print("Tangentials (no-slip ~0; free-slip unconstrained):")
    print(f"  max |v_tan| on E coast : {max_vt_E:0.3e}  (rel={rel_vt_E:0.3e})")
    print(f"  max |u_tan| on N coast : {max_ut_N:0.3e}  (rel={rel_ut_N:0.3e})")
    print("Normal derivative of tangential (free-slip ~0; index units):")
    print(f"  max |∂v/∂n| on E coast : {max_dv_dn:0.3e} (index units)")
    print(f"  max |∂u/∂n| on N coast : {max_du_dn:0.3e} (index units)")

    # Heuristic PASS/FAIL (soft) with NaN-safe logic
    # Tolerances: normals must be < 1e-3 of interior scale, tangentials for noslip also < 1e-3,
    # derivatives for free/noratio must be < 1e-3 of interior scale in index units.
    # If scales are zero, skip PASS/FAIL (cannot judge).
    tol_rel = 1e-3
    def rel_ok(val, scale):
        if not np.isfinite(val):
            return False
        if scale <= 0 or not np.isfinite(scale):
            return True  # cannot judge, don't fail
        return (val <= tol_rel * scale)

    normals_ok = rel_ok(max_un_E, max(vE_ref, uN_ref)) and rel_ok(max_vn_N, max(vE_ref, uN_ref))
    noslip_ok  = rel_ok(max_vt_E, vE_ref) and rel_ok(max_ut_N, uN_ref)
    freeslip_ok = rel_ok(max_dv_dn, vE_ref) and rel_ok(max_du_dn, uN_ref)

    if expect == "noslip":
        print(f"Verdict (no-slip): normals={'OK' if normals_ok else 'FAIL'}, "
              f"tangentials={'OK' if noslip_ok else 'FAIL'}")
    elif expect in ("free", "noratio"):
        print(f"Verdict (free-slip): normals={'OK' if normals_ok else 'FAIL'}, "
              f"d(tangent)/dn={'OK' if freeslip_ok else 'FAIL'}  "
              f"(note: tangentials may be O(1) of interior and are NOT a failure)")
    else:
        print(f"Verdict: normals={'OK' if normals_ok else 'FAIL'}; "
              f"noslip tangentials={'OK' if noslip_ok else 'FAIL'}; "
              f"free-slip d(tan)/dn={'OK' if freeslip_ok else 'FAIL'}")

    print("")  # blank line
    ds.close()


def main():
    ap = argparse.ArgumentParser(description="Check CICE BCs (impermeability, no-slip, free-slip) from instantaneous files.")
    ap.add_argument("--exp", action="append", default=[],
                    help="Experiment in the form LABEL:PATH where LABEL in {free,noratio,noslip,any}. "
                         "PATH is a history/ directory containing iceh_inst.*.nc. Repeatable.")
    ap.add_argument("--file", default=None,
                    help="Explicit file name to analyze inside each PATH; if omitted, newest iceh_inst.*.nc is used.")
    args = ap.parse_args()

    if not args.exp:
        print("Provide at least one --exp LABEL:PATH. Example:")
        print("  --exp free:/g/data/.../BC_free/history --exp noratio:/g/data/.../BC_free_noratio/history --exp noslip:/g/data/.../BC_noslip/history")
        sys.exit(2)

    exps: Dict[str, Tuple[str, str]] = {}  # label -> (path, expect)
    for item in args.exp:
        if ":" not in item:
            print(f"Bad --exp '{item}'. Expected LABEL:PATH")
            sys.exit(2)
        label, path = item.split(":", 1)
        label = label.strip()
        path = path.strip()
        if not os.path.isdir(path):
            print(f"[{label}] Not a directory: {path}")
            continue
        # Map label to expected behavior
        expect = "noslip" if label.lower() in ("noslip", "no-slip") else \
                 "free"   if label.lower() in ("free", "free-slip") else \
                 "noratio" if label.lower() in ("noratio", "free-slip-noratio") else "any"
        exps[label] = (path, expect)

    if not exps:
        print("No valid experiments found.")
        sys.exit(1)

    for label, (path, expect) in exps.items():
        if args.file:
            nc_path = os.path.join(path, args.file)
            if not os.path.isfile(nc_path):
                print(f"[{label}] File not found: {nc_path}")
                continue
        else:
            nc_path = find_latest_inst_file(path)
            if nc_path is None:
                print(f"[{label}] No iceh_inst.*.nc found under {path}")
                continue
        analyze_one_file(nc_path, label, expect)

if __name__ == "__main__":
    main()
