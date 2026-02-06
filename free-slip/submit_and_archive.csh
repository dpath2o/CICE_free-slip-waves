#!/bin/csh -f
# ------------------------------------------------------------
# submit_ana_archive.csh
#
# Runs CICE in <=2-year chunks (as per your existing ice_in stop settings),
# auto-detects initial vs restart, enforces restart consistency, and chains
# sequential qsub runs until END_DATE achieved or a crash is detected.
#
# Usage:
#   ./submit_ana_archive.csh -s SIM_NAME -S YYYY-MM-DD -E YYYY-MM-DD [-c CASE_NAME]
#
# Notes:
# - If -c omitted, CASE_NAME = basename($PWD)
# - Requires: python3, perl, qsub, qstat on PATH.
# ------------------------------------------------------------

set nonomatch

# -------- args --------
set CASE  = ""
set SIM   = ""
set START = ""
set END   = ""

if ( $#argv == 0 ) goto usage

while ( $#argv > 0 )
  switch ( "$1" )
    case "-c":
    case "--case":
      shift; if ( $#argv == 0 ) goto usage
      set CASE = "$1"
      breaksw
    case "-s":
    case "--sim":
      shift; if ( $#argv == 0 ) goto usage
      set SIM = "$1"
      breaksw
    case "-S":
    case "--start":
      shift; if ( $#argv == 0 ) goto usage
      set START = "$1"
      breaksw
    case "-E":
    case "--end":
      shift; if ( $#argv == 0 ) goto usage
      set END = "$1"
      breaksw
    case "-h":
    case "--help":
      goto usage
    default:
      echo "ERROR: Unknown argument: $1"
      goto usage
  endsw
  shift
end

if ( "$CASE" == "" ) set CASE = `basename "$PWD"`
if ( "$SIM"  == "" || "$START" == "" || "$END" == "" ) then
  echo "ERROR: Missing required args."
  goto usage
endif

# -------- constants / paths --------
set CHUNK_YEARS = 2
set SLEEP_S     = 60

set RUNROOT     = "$HOME/CICE_runs/$CASE"
set ARCH        = "$HOME/AFIM_archive/$SIM"
set ARCH_HIST   = "$ARCH/history/daily"
set ARCH_RST    = "$ARCH/restart"
set FSLIPROOT   = "$HOME/CICE_free-slip/$CASE"

# wrappers
set CLASS_DIR   = "$HOME/AFIM/src/AFIM/scripts/classification"
set METRICS_DIR = "$HOME/AFIM/src/AFIM/scripts/metrics"
set METLOG      = "$HOME/logs/metrics_${SIM}_ispd_thresh0.0005.log"

# -------- sanity --------
if ( ! -d "$RUNROOT" ) then
  echo "ERROR: Missing run directory: $RUNROOT"
  exit 2
endif
if ( ! -e "$FSLIPROOT/ice_in" ) then
  echo "ERROR: Missing $FSLIPROOT/ice_in"
  exit 2
endif
if ( ! -e "$FSLIPROOT/cice.run" ) then
  echo "ERROR: Missing $FSLIPROOT/cice.run"
  exit 2
endif
if ( ! -d "$RUNROOT/restart" ) then
  echo "ERROR: Missing restart dir: $RUNROOT/restart"
  exit 2
endif
if ( ! -d "$RUNROOT/history" ) then
  echo "ERROR: Missing history dir: $RUNROOT/history"
  exit 2
endif

mkdir -p "$ARCH_HIST" "$ARCH_RST" "$ARCH"

# -------- helper: set ice_in flags safely --------
# (uses perl regex; assumes keys exist in ice_in)
alias set_runtype_continue  'perl -pi -e "s/^(\\s*runtype\\s*=\\s*).*/\\1\\x27continue\\x27/i" "$RUNROOT/ice_in"; perl -pi -e "s/^(\\s*use_restart_time\\s*=\\s*).*/\\1.true./i" "$RUNROOT/ice_in"'
alias set_runtype_initial   'perl -pi -e "s/^(\\s*runtype\\s*=\\s*).*/\\1\\x27initial\\x27/i" "$RUNROOT/ice_in";  perl -pi -e "s/^(\\s*use_restart_time\\s*=\\s*).*/\\1.false./i" "$RUNROOT/ice_in"'

# stop controls: prefer nyears when aligned; otherwise ndays
alias set_stop_nyears 'perl -pi -e "s/^(\\s*stop_option\\s*=\\s*).*/\\1\\x27nyears\\x27/i" "$RUNROOT/ice_in"; perl -pi -e "s/^(\\s*stop_n\\s*=\\s*).*/\\1$STOP_N/i" "$RUNROOT/ice_in"'
alias set_stop_ndays  'perl -pi -e "s/^(\\s*stop_option\\s*=\\s*).*/\\1\\x27ndays\\x27/i"  "$RUNROOT/ice_in"; perl -pi -e "s/^(\\s*stop_n\\s*=\\s*).*/\\1$STOP_N/i" "$RUNROOT/ice_in"'

# -------- helper: parse date from iced.YYYY-MM-DD-00000.nc --------
set restart_ptr = ""
set restart_date_ptr = ""

if ( -e "$RUNROOT/ice.restart_file" ) then
  set restart_ptr = `head -1 "$RUNROOT/ice.restart_file" | tr -d ' \t\r\n'`
  if ( "$restart_ptr" != "" && -e "$RUNROOT/restart/$restart_ptr" ) then
    set restart_date_ptr = `python3 - << 'PY'
import re,sys
fn = sys.argv[1]
m = re.search(r"iced\.(\d{4}-\d{2}-\d{2})-00000\.nc", fn)
print(m.group(1) if m else "")
PY "$restart_ptr"`
  endif
endif

# -------- decide initial vs restart --------
set MODE = "initial"
set CUR_START = "$START"

if ( "$restart_date_ptr" != "" ) then
  # restart is valid and consistent
  set MODE = "restart"
  set CUR_START = "$restart_date_ptr"
endif

echo "`date` submit_ana_archive.csh: CASE=$CASE SIM=$SIM MODE=$MODE START=$START END=$END" >> "$RUNROOT/README.case"
echo "`date` MODE=$MODE; starting chunk loop from $CUR_START"

# If restarting, enforce your hard rules:
if ( "$MODE" == "restart" ) then
  if ( ! -e "$RUNROOT/ice.restart_file" ) then
    echo "ERROR: restart mode but missing $RUNROOT/ice.restart_file"
    exit 3
  endif
  if ( "$restart_ptr" == "" ) then
    echo "ERROR: restart mode but ice.restart_file is empty"
    exit 3
  endif
  if ( ! -e "$RUNROOT/restart/$restart_ptr" ) then
    echo "ERROR: restart mode but missing $RUNROOT/restart/$restart_ptr"
    exit 3
  endif
  set_runtype_continue
endif

# If initialising, enforce initial flags for first run only:
if ( "$MODE" == "initial" ) then
  set_runtype_initial
endif

# -------- main loop --------
set LAST_GOOD_END = ""

while ( "$CUR_START" <= "$END" )

  # compute chunk end, restart date, and STOP_N / STOP_OPTION choice
  # Outputs: seg_end restart_date stop_opt stop_n
  set py = (`python3 - << 'PY'
from datetime import date, timedelta
import sys

def parse(s):
    y,m,d = map(int, s.split("-"))
    return date(y,m,d)

def add_years(d, years):
    try:
        return d.replace(year=d.year + years)
    except ValueError:
        # handle Feb 29 -> Feb 28
        return d.replace(month=2, day=28, year=d.year + years)

cur = parse(sys.argv[1])
end = parse(sys.argv[2])
chunk_years = int(sys.argv[3])

# nominal chunk end = (cur + chunk_years years) - 1 day
nom_end = add_years(cur, chunk_years) - timedelta(days=1)
seg_end = nom_end if nom_end <= end else end
restart_date = seg_end + timedelta(days=1)

# choose stop control:
# if cur is Jan-01 and restart_date is Jan-01 => nyears is safe
stop_opt = "ndays"
stop_n = (restart_date - cur).days  # number of days to advance to restart_date

if cur.month == 1 and cur.day == 1 and restart_date.month == 1 and restart_date.day == 1:
    years = restart_date.year - cur.year
    if years >= 1:
        stop_opt = "nyears"
        stop_n = years

print(seg_end.isoformat(), restart_date.isoformat(), stop_opt, stop_n)
PY "$CUR_START" "$END" "$CHUNK_YEARS"`)

  set SEG_END      = "$py[1]"
  set RESTART_DATE = "$py[2]"
  set STOP_OPT     = "$py[3]"
  set STOP_N       = "$py[4]"

  echo "`date` chunk: CUR_START=$CUR_START SEG_END=$SEG_END RESTART_DATE=$RESTART_DATE STOP_OPT=$STOP_OPT STOP_N=$STOP_N"

  # set stop controls in ice_in for this chunk
  if ( "$STOP_OPT" == "nyears" ) then
    set_stop_nyears
  else
    set_stop_ndays
  endif

  # Ensure correct runtype flags for this chunk
  if ( "$MODE" == "initial" ) then
    # first chunk only
    set_runtype_initial
  else
    # restart/continue chunks
    # enforce restart pointer validity before submitting
    if ( ! -e "$RUNROOT/ice.restart_file" ) then
      echo "ERROR: continue mode but missing $RUNROOT/ice.restart_file"
      exit 4
    endif
    set restart_ptr = `head -1 "$RUNROOT/ice.restart_file" | tr -d ' \t\r\n'`
    if ( "$restart_ptr" == "" ) then
      echo "ERROR: continue mode but ice.restart_file empty"
      exit 4
    endif
    if ( ! -e "$RUNROOT/restart/$restart_ptr" ) then
      echo "ERROR: continue mode but missing $RUNROOT/restart/$restart_ptr"
      exit 4
    endif
    set_runtype_continue
  endif

  # submit run
  cd "$RUNROOT" || exit 5
  set jobid = `qsub ./cice.run`
  if ( "$jobid" == "" ) then
    echo "ERROR: qsub returned empty job id"
    exit 6
  endif
  echo "`date` submit_ana_archive.csh: submitted $jobid for chunk $CUR_START -> $SEG_END" >> "$RUNROOT/README.case"
  echo "`date` waiting for PBS job $jobid..."

  # wait for completion
  while ( 1 )
    qstat "$jobid" >& /dev/null
    if ( $status != 0 ) break
    sleep $SLEEP_S
  end

  # detect success: expected restart file must exist
  set EXPECTED_RST = "iced.${RESTART_DATE}-00000.nc"
  if ( ! -e "$RUNROOT/restart/$EXPECTED_RST" ) then
    echo "`date` ERROR: expected restart missing: $RUNROOT/restart/$EXPECTED_RST"
    echo "`date` Treating as MODEL CRASH; stopping chain."

    # archive what we can (logs/diag) to inspect
    if ( -e "$RUNROOT/ice_diag.d" ) mv "$RUNROOT/ice_diag.d" "$ARCH"/
    set runlogs = ( "$RUNROOT"/cice.runlog.* )
    if ( "$runlogs[1]" != "$RUNROOT/cice.runlog.*" ) mv $runlogs "$ARCH"/

    # if you want an explicit crash marker:
    echo "`date` CRASH during chunk starting $CUR_START (expected $EXPECTED_RST)" >> "$ARCH/CRASH_MARKER.txt"

    # run analysis up to last good end (if any)
    if ( "$LAST_GOOD_END" != "" ) then
      echo "`date` Running AFIM wrappers for achieved period: $START -> $LAST_GOOD_END"
      cd "$CLASS_DIR" || exit 7
      ./classify_sea_ice_pbs_wrapper.sh -s "$SIM" -i FI -g Tc -S "$START" -E "$LAST_GOOD_END" -z -y -k
      cd "$METRICS_DIR" || exit 7
      ./metrics_pbs_wrapper.sh -s "$SIM" -i FI -g Tc -S "$START" -E "$LAST_GOOD_END" -z -p
      if ( -e "$METLOG" ) tail "$METLOG"
    endif

    exit 9
  endif

  # optional “fatal” scan (belt-and-braces)
  if ( -e "$RUNROOT/ice_diag.d" ) then
    grep -Eiq "FATAL|ABORT|SIGSEGV|MPI_ABORT|ERROR STOP" "$RUNROOT/ice_diag.d"
    if ( $status == 0 ) then
      echo "`date` WARNING: fatal markers found in ice_diag.d; stopping chain (restart exists but logs look bad)."
      mv "$RUNROOT/ice_diag.d" "$ARCH"/
      exit 10
    endif
  endif

  # mark successful chunk
  set LAST_GOOD_END = "$SEG_END"

  # -------- archive step (per chunk) --------
  mkdir -p "$ARCH_HIST" "$ARCH_RST" "$ARCH"

  # move history
  set histfiles = ( "$RUNROOT/history"/iceh.* )
  if ( "$histfiles[1]" != "$RUNROOT/history/iceh.*" ) then
    mv $histfiles "$ARCH_HIST"/
  endif

  # move run logs
  set runlogs = ( "$RUNROOT"/cice.runlog.* )
  if ( "$runlogs[1]" != "$RUNROOT/cice.runlog.*" ) then
    mv $runlogs "$ARCH"/
  endif

  # move ice_diag.d (keep ice_in in place for next chunk!)
  if ( -e "$RUNROOT/ice_diag.d" ) mv "$RUNROOT/ice_diag.d" "$ARCH"/

  # move free-slip logs
  set fsliplogs = ( "$FSLIPROOT"/free-slip.o* )
  if ( "$fsliplogs[1]" != "$FSLIPROOT/free-slip.o*" ) then
    mv $fsliplogs "$ARCH"/
  endif

  # move *the* restart file to archive, but keep a symlink in run directory
  if ( -e "$RUNROOT/restart/$EXPECTED_RST" ) then
    # remove any stale archived copy
    if ( -e "$ARCH_RST/$EXPECTED_RST" ) rm -f "$ARCH_RST/$EXPECTED_RST"
    mv "$RUNROOT/restart/$EXPECTED_RST" "$ARCH_RST"/
    ln -sf "$ARCH_RST/$EXPECTED_RST" "$RUNROOT/restart/$EXPECTED_RST"
  endif

  # remove other restart files (keep the expected one)
  foreach f ( "$RUNROOT/restart"/iced.* )
    if ( "$f" != "$RUNROOT/restart/iced.*" ) then
      if ( "`basename "$f"`" != "$EXPECTED_RST" ) rm -f "$f"
    endif
  end

  # update ice.restart_file pointer for next run
  echo "$EXPECTED_RST" > "$RUNROOT/ice.restart_file"

  # after first successful run, switch mode to restart for subsequent chunks
  set MODE = "restart"
  set CUR_START = "$RESTART_DATE"

end

# -------- reached END successfully --------
echo "`date` Completed all chunks through END=$END (last_good_end=$LAST_GOOD_END)"
echo "`date` Running AFIM wrappers for full period: $START -> $END"

cd "$CLASS_DIR" || exit 11
./classify_sea_ice_pbs_wrapper.sh -s "$SIM" -i FI -g Tc -S "$START" -E "$END" -z -y -k

cd "$METRICS_DIR" || exit 11
./metrics_pbs_wrapper.sh -s "$SIM" -i FI -g Tc -S "$START" -E "$END" -z -p

if ( -e "$METLOG" ) tail "$METLOG"
exit 0

usage:
echo "Usage:"
echo "  $0 -s SIM_NAME -S YYYY-MM-DD -E YYYY-MM-DD [-c CASE_NAME]"
echo ""
echo "Behaviour:"
echo "  - Runs in <=2-year chunks until END or crash."
echo "  - Auto-detects restart if ice.restart_file + matching restart/iced.* present."
echo "  - Enforces: runtype/use_restart_time consistent with initial vs continue."
exit 1