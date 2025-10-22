#!/usr/bin/env bash
set -euo pipefail

# nc_equiv_check.sh  A.nc  B.nc  [--abs-tol 1e-12] [--rel-tol 1e-9] [--no-edges]
# Compare two NetCDF files for numerical equivalence using NCO only.

ABS_TOL="1e-12"
REL_TOL="1e-9"
DO_EDGES=1

# Variables to consider (both with and without _1)
VARS=(uvelE_1 vvelE_1 uvelN_1 vvelN_1 uvel_1 vvel_1 divu_1 vort_1
      uvelE   vvelE   uvelN   vvelN   uvel   vvel   divu   vort  )

usage(){ echo "Usage: $0 A.nc B.nc [--abs-tol X] [--rel-tol Y] [--no-edges]"; exit 2; }

[[ $# -lt 2 ]] && usage
A_IN=$1; shift
B_IN=$1; shift

while [[ $# -gt 0 ]]; do
  case "$1" in
    --abs-tol) ABS_TOL="$2"; shift 2;;
    --rel-tol) REL_TOL="$2"; shift 2;;
    --no-edges) DO_EDGES=0; shift;;
    *) echo "Unknown arg: $1"; usage;;
  esac
done

# first check bitwise
cmp -s "$A_IN" "$B_IN" && echo "BITWISE IDENTICAL" || echo "bytes differ"

# perform MDF5 sum
 md5sum "$A_IN" "$B_IN"

need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: need $1"; exit 1; }; }
for exe in ncks ncdiff ncap2 ncwa ncatted ncrename; do need "$exe"; done

tmp=$(mktemp -d -t nc_equiv.XXXXXX)
trap 'rm -rf "$tmp"' EXIT

A_TMP="$tmp/A.tmp.nc"; B_TMP="$tmp/B.tmp.nc"
A_STD="$tmp/A.std.nc"; B_STD="$tmp/B.std.nc"
DIFF="$tmp/diff.nc"; DIFF_ABS="$tmp/diff_abs.nc"; DIFF_MAX="$tmp/diff_maxabs.nc"
A_MAX="$tmp/A_max.nc"; REL="$tmp/rel_max.nc"
DIMS="$tmp/dims.nc"

cp -f "$A_IN" "$A_TMP"
cp -f "$B_IN" "$B_TMP"

# Strip 'history' (ignore failures if attribute absent)
ncatted -O -a history,global,d,, "$A_TMP" || true
ncatted -O -a history,global,d,, "$B_TMP" || true

# Normalise to classic, no compression; fix record dim if present
normalise() {
  local IN="$1" OUT="$2"
  if ncks -O -3 -L 0 --fix_rec_dmn time "$IN" "$OUT" 2>/dev/null; then
    :
  else
    ncks -O -3 -L 0 "$IN" "$OUT"
  fi
}
normalise "$A_TMP" "$A_STD"
normalise "$B_TMP" "$B_STD"

# bitwise after NC normalisation 
cmp -s "$A_STD" "$B_STD" && echo "IDENTICAL AFTER NORMALISE" || echo "still different"

echo "=== Header diff (metadata/dims/vars) ==="
diff -u <(ncks -m "$A_STD") <(ncks -m "$B_STD") || true
echo "========================================"

# Numeric difference
ncdiff -O "$B_STD" "$A_STD" "$DIFF"

# Build list of variables actually present in DIFF by probing ncks exit status
present_vars=()
for v in "${VARS[@]}"; do
  if ncks -q -C -H -v "$v" "$DIFF" >/dev/null 2>&1; then
    present_vars+=("$v")
  fi
done

if [[ ${#present_vars[@]} -eq 0 ]]; then
  echo "No target variables present to compare. Exiting."
  exit 0
fi

# abs() of differences (NaNs remain NaN)
ABS_EXPR=""
for v in "${present_vars[@]}"; do ABS_EXPR+="$v=abs($v);"; done
ncap2 -O -s "$ABS_EXPR" "$DIFF" "$DIFF_ABS"

# Reduce to global maxima over time,nj,ni
ncwa -O -y max -a time,nj,ni "$DIFF_ABS" "$DIFF_MAX"

# Reference maxima from A with the same variables
ncwa -O -y max -a time,nj,ni -v "$(IFS=,; echo "${present_vars[*]}")" "$A_STD" "$A_MAX"
for v in "${present_vars[@]}"; do ncrename -O -v "$v","${v}_ref" "$A_MAX"; done
ncks -A "$A_MAX" "$DIFF_MAX"

# Relative ratios
REL_EXPR="eps=${ABS_TOL};"
for v in "${present_vars[@]}"; do REL_EXPR+="${v}_rel=${v}/(${v}_ref+eps);"; done
ncap2 -O -s "$REL_EXPR" "$DIFF_MAX" "$REL"

echo
echo "=== Numeric comparison (max over time,nj,ni) ==="
printf "%-12s %16s %16s %16s %8s\n" "variable" "max|delta|" "max|ref|" "rel=delta/(ref+eps)" "PASS?"
all_pass=1
for v in "${present_vars[@]}"; do
  d=$(ncks -H -C -v "$v"        "$REL" | awk -F'= ' '/=/{print $2}' | tr -d ' ;')
  r=$(ncks -H -C -v "${v}_ref"  "$REL" | awk -F'= ' '/=/{print $2}' | tr -d ' ;')
  q=$(ncks -H -C -v "${v}_rel"  "$REL" | awk -F'= ' '/=/{print $2}' | tr -d ' ;')
  pass=$(awk -v d="${d:-nan}" -v q="${q:-nan}" -v at="$ABS_TOL" -v rt="$REL_TOL" 'BEGIN{
    ok=0;
    if (d==d && d<=at) ok=1;           # finite and <= abs tol
    if (q==q && q<=rt) ok=1;           # finite and <= rel tol
    print ok ? "OK" : "FAIL";
  }')
  printf "%-12s %16.6e %16.6e %16.6e %8s\n" "$v" "${d:-nan}" "${r:-nan}" "${q:-nan}" "$pass"
  [[ "$pass" == "OK" ]] || all_pass=0
done

# -------- Edge sanity check (coastal normals only) --------
echo
echo "=== Edge sanity check (coastal normals only) ==="

# Pick available face-normal names (either with or without _1), from the diff file
pick_var() {
  local a="$1" b="$2"
  if ncks -q -C -H -v "$a" "$DIFF_ABS" >/dev/null 2>&1; then echo "$a"
  elif ncks -q -C -H -v "$b" "$DIFF_ABS" >/dev/null 2>&1; then echo "$b"
  else echo ""
  fi
}

E_VAR=$(pick_var uvelE_1 uvelE)   # normal on E-faces
N_VAR=$(pick_var vvelN_1 vvelN)   # normal on N-faces

# Robust ni/nj from metadata (works with ncks everywhere)
get_dim_size() {
  local file="$1" dim="$2"
  ncks -m "$file" \
  | awk -v d="$dim" '
      BEGIN{ IGNORECASE=1 }
      $1==d && $2=="=" { gsub(/;/,"",$3); print $3; exit }
    '
}

NI=$(get_dim_size "$A_STD" ni || true)
NJ=$(get_dim_size "$A_STD" nj || true)
NI=$(($NI-1)) # index adjustment for zeroth counting
NJ=$(($NJ-1))

if [[ -z "${NI:-}" || -z "${NJ:-}" ]]; then
  echo "Could not determine ni/nj from $A_STD; skipping edge check."
else
  # west/east edges for E faces (max over time,nj)
  if [[ -n "$E_VAR" ]]; then
    ncwa -O -y max -a time,nj -d ni,1,1          -v "$E_VAR" "$DIFF_ABS" "$tmp/max_uE_west.nc"
    ncwa -O -y max -a time,nj -d ni,"$NI","$NI"  -v "$E_VAR" "$DIFF_ABS" "$tmp/max_uE_east.nc"
    max_w=$(ncks -H -C -v "$E_VAR" "$tmp/max_uE_west.nc" | awk -F'= ' '/=/{gsub(/[ ;]/,"",$2); print $2}')
    max_e=$(ncks -H -C -v "$E_VAR" "$tmp/max_uE_east.nc" | awk -F'= ' '/=/{gsub(/[ ;]/,"",$2); print $2}')
    echo "max |delta_uN| on E faces: west=${max_w:-nan}, east=${max_e:-nan}"
  else
    echo "E-face normal variable (uvelE_1/uvelE) not found in diff; skipping E edge check."
  fi

  # south/north edges for N faces (max over time,ni)
  if [[ -n "$N_VAR" ]]; then
    ncwa -O -y max -a time,ni -d nj,1,1          -v "$N_VAR" "$DIFF_ABS" "$tmp/max_vN_south.nc"
    ncwa -O -y max -a time,ni -d nj,"$NJ","$NJ"  -v "$N_VAR" "$DIFF_ABS" "$tmp/max_vN_north.nc"
    max_s=$(ncks -H -C -v "$N_VAR" "$tmp/max_vN_south.nc" | awk -F'= ' '/=/{gsub(/[ ;]/,"",$2); print $2}')
    max_n=$(ncks -H -C -v "$N_VAR" "$tmp/max_vN_north.nc" | awk -F'= ' '/=/{gsub(/[ ;]/,"",$2); print $2}')
    echo "max |delta_vN| on N faces: south=${max_s:-nan}, north=${max_n:-nan}"
  else
    echo "N-face normal variable (vvelN_1/vvelN) not found in diff; skipping N edge check."
  fi
fi

echo
if [[ $all_pass -eq 1 ]]; then
  echo "OVERALL: PASS (within |delta| <= ${ABS_TOL} OR relative <= ${REL_TOL})"
  exit 0
else
  echo "OVERALL: FAIL (at least one variable exceeds tolerances)"
  exit 1
fi
