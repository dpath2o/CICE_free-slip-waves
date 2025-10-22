#!/bin/csh -f

# must load an appropriate python environment;
# I'm using xp65 on gadi because it's stable and good
if (! $?MODULESHOME) then
  source /etc/profile.d/modules.csh
endif
module use /g/data/xp65/public/modules
module load conda/analysis3-25.05

./CDP_analysis.py

# setup
set case = "BCtest1"
set casedir = "/g/data/gv90/da1339/cice-dirs/runs/BCtest1/"
set histdir = "${casedir}/history"
set files = ("${histdir}/iceh.2005-01-01.nc" \
             "${histdir}/iceh.2005-01-02.nc" \
             "${histdir}/iceh.2005-01-03.nc" \
             "${histdir}/iceh.2005-01-04.nc" \
             "${histdir}/iceh.2005-01-05.nc" )
set notes = ("Coastal Drag test: 01 Jan 05" \
             "Coastal Drag test: 02 Jan 05" \
             "Coastal Drag test: 03 Jan 05" \
             "Coastal Drag test: 04 Jan 05" \
             "Coastal Drag test: 05 Jan 05" )
set fstrs = ("01Jan05" \
             "02Jan05" \
             "03Jan05" \
             "04Jan05" \
             "05Jan05" )
set fields = ("aice" "hi" "KuxN" "KuyN" "KuxE" "KuyE" "uvel" "vvel")
mkdir -pv figs
# call the python plotting routines
echo " "
echo " "
echo ./timeseries.py \"${casedir}\" --case \"${case}\" --grid
./timeseries.py "${casedir}" --case "${case}" --grid
echo " "
set cnt = 0
while ($cnt < ${#files})
  @ cnt = $cnt + 1
  set file = "${files[$cnt]}"
  set note = "${notes[$cnt]}"
  set fstr = "${fstrs[$cnt]}"
  foreach field ($fields)
    echo ./ciceplots2d.py \"$field\" \"$file\" \"$case\" \"$note\" \"$fstr\"
    ./ciceplots2d.py "$field" "$file" "$case" "$note" "$fstr"
  end
end
echo "DONE"

