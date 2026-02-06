#!/bin/csh -f
# Usage:
#   ./submit_and_archive.csh -s SIM_NAME -S YYYY-MM-DD -E YYYY-MM-DD [-c CASE_NAME]
# If -c not given: CASE_NAME = basename($PWD)

set CASE_NAME   = ""
set SIM_NAME    = ""
set START_DATE  = ""
set END_DATE    = ""

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
      echo "Usage: $0 -s SIM_NAME -S YYYY-MM-DD -E YYYY-MM-DD [-c CASE_NAME]"
      exit 2
  endsw
  shift
end

if ( "$SIM_NAME" == "" || "$START_DATE" == "" || "$END_DATE" == "" ) then
  echo "ERROR: must provide -s SIM_NAME -S START_DATE -E END_DATE"
  echo "Usage: $0 -s SIM_NAME -S YYYY-MM-DD -E YYYY-MM-DD [-c CASE_NAME]"
  exit 2
endif

if ( "$CASE_NAME" == "" ) then
  set CASE_NAME = `basename "$PWD"`
endif

# Where cice.run will write its chain status flags
set ARCHROOT  = "$HOME/AFIM_archive/$SIM_NAME"
mkdir -p "$ARCHROOT"

set DONE_FLAG  = "$ARCHROOT/cice_chain_DONE_${START_DATE}_${END_DATE}.flag"
set CRASH_FLAG = "$ARCHROOT/cice_chain_CRASH_${START_DATE}_${END_DATE}.flag"

# Export vars so cice.run knows what it is aiming for, and where to signal status
setenv CASE_NAME   "$CASE_NAME"
setenv SIM_NAME    "$SIM_NAME"
setenv START_DATE  "$START_DATE"
setenv END_DATE    "$END_DATE"
setenv DONE_FLAG   "$DONE_FLAG"
setenv CRASH_FLAG  "$CRASH_FLAG"

# Submit ONE job only. cice.run is responsible for chaining 2-year chunks + restart logic.
set first_jobid = "`qsub -V ./cice.run`"
if ( "$first_jobid" == "" ) then
  echo "ERROR: qsub returned empty job id."
  exit 2
endif

echo "`date` submit_and_archive: submitted CICE chain start job=$first_jobid CASE=$CASE_NAME SIM=$SIM_NAME ($START_DATE..$END_DATE)" > ${PWD}/README.case
echo "Waiting for DONE flag: $DONE_FLAG" >> ${PWD}/README.case
echo "Crash flag is:          $CRASH_FLAG" >> ${PWD}/README.case

# Wait for CICE chain to finish (DONE) or fail (CRASH)
while ( ! -e "$DONE_FLAG" )
  if ( -e "$CRASH_FLAG" ) then
    echo "ERROR: CICE chain crashed."
    echo "See: $CRASH_FLAG"
    exit 1
  endif
  sleep 600
end

echo "`date` submit_and_archive: CICE chain finished OK: `cat $DONE_FLAG`" >> ${PWD}/README.case

# -----------------------------
# Post-processing (only after DONE):
#   1) fast-ice classification (PBS jobs)
#   2) metrics (PBS) AFTER classification completes successfully
# -----------------------------

set CLASS_WRAP = "$HOME/AFIM/src/AFIM/scripts/classification/classify_sea_ice_pbs_wrapper.sh"
set MET_PBS    = "$HOME/AFIM/src/AFIM/scripts/metrics/metrics.pbs"

# Submit classification (yearly) and capture job IDs robustly
set CLASS_OUT = "$ARCHROOT/classify_submit_${START_DATE}_${END_DATE}_$$.out"
echo "`date` submit_and_archive: submitting classification..." >> ${PWD}/README.case

bash "$CLASS_WRAP" -s "$SIM_NAME" -i FI -g Tc -S "$START_DATE" -E "$END_DATE" -z -y -k |& tee "$CLASS_OUT"

set class_jobids = ()
foreach jid ( `awk '/^[0-9]/{print $1}' "$CLASS_OUT"` )
  set class_jobids = ( $class_jobids $jid )
end

if ( $#class_jobids == 0 ) then
  echo "ERROR: No classification job IDs captured. See: $CLASS_OUT"
  exit 2
endif

set dep = `echo $class_jobids | tr ' ' ':'`
echo "Classification job IDs: $class_jobids" >> ${PWD}/README.case
echo "Metrics dependency: afterok:$dep" >> ${PWD}/README.case

# Build env file for metrics.pbs (ENVFILE contract)
set MET_ENV = "$HOME/logs/ENVFILE_metrics_${SIM_NAME}_$$.sh"
cat >! "$MET_ENV" << EOF
#!/bin/bash
export sim_name="$SIM_NAME"
export ispd_thresh="0.0005"
export ice_type="FI"
export borc2t_type="Tc"
export start_date="$START_DATE"
export end_date="$END_DATE"
export overwrite_zarr=true
export overwrite_png=true
export compute_boolean=false
EOF
chmod +x "$MET_ENV"

# Submit metrics job AFTER ALL classification jobs succeed
set met_jobid = "`qsub -W depend=afterok:$dep -v ENVFILE="$MET_ENV" "$MET_PBS"`"
if ( "$met_jobid" == "" ) then
  echo "ERROR: metrics qsub returned empty job id."
  exit 2
endif

echo "Submitted metrics job: $met_jobid" >> ${PWD}/README.case
echo "Tail metrics log:" >> ${PWD}/README.case
echo "  tail $HOME/logs/metrics_${SIM_NAME}_ispd_thresh0.0005.log" >> ${PWD}/README.case

echo "`date` submit_and_archive: DONE (submitted post-processing)" >> ${PWD}/README.case