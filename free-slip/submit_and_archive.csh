#!/bin/csh -f
# Usage:
#   ./submit_and_archive.csh -s SIM_NAME -S YYYY-MM-DD -E YYYY-MM-DD [-c CASE_NAME]
#
# If -c not given: CASE_NAME = basename($PWD)

set CASE_NAME = ""
set SIM_NAME  = ""
set START_DATE = ""
set END_DATE   = ""

while ( $#argv )
  switch ( "$argv[1]" )
    case "-c":
      shift
      set CASE_NAME = "$argv[1]"
      breaksw
    case "-s":
      shift
      set SIM_NAME = "$argv[1]"
      breaksw
    case "-S":
      shift
      set START_DATE = "$argv[1]"
      breaksw
    case "-E":
      shift
      set END_DATE = "$argv[1]"
      breaksw
    default:
      echo "Unknown arg: $argv[1]"
      exit 2
  endsw
  shift
end

if ( "$SIM_NAME" == "" || "$START_DATE" == "" || "$END_DATE" == "" ) then
  echo "ERROR: must provide -s SIM_NAME -S START_DATE -E END_DATE"
  exit 2
endif

if ( "$CASE_NAME" == "" ) then
  set CASE_NAME = `basename "$PWD"`
endif

# Export vars for the PBS job
setenv CASE_NAME   "$CASE_NAME"
setenv SIM_NAME    "$SIM_NAME"
setenv START_DATE  "$START_DATE"
setenv END_DATE    "$END_DATE"

# Submit the first chunk. cice.run will auto-chain subsequent chunks.
qsub -V ./cice.run

echo "`date` ${0}: submitted first chunk for CASE=$CASE_NAME SIM=$SIM_NAME ($START_DATE .. $END_DATE)" >> ${PWD}/README.case

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