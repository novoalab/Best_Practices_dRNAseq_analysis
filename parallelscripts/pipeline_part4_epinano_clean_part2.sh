#!/bin/bash


INPUT=$1

source /users/enovoa/joramirez/.bashrc
 
for d in ${INPUT}/ALBACORE_2.1.7_OUTPUT ${INPUT}/ALBACORE_2.3.4_OUTPUT ${INPUT}/GUPPY_2.3.1_OUTPUT ${INPUT}/GUPPY_3.0.3_OUTPUT; do
if [[ ! -f $d/epinano/minimap2/logs2/qout ]]
then
 mkdir -p $d/epinano/minimap2/logs2/qout
 mkdir -p $d/epinano/minimap2/logs2/qerr
 mkdir -p $d/epinano/graphmap/logs2/qout
 mkdir -p $d/epinano/graphmap/logs2/qerr
else
 echo "[ WARNING ] Found existing files. Clean up?"
fi
NUMJOBS_m=$( ls $d/epinano/minimap2/tsv/ | wc -l )  
NUMJOBS_g=$( ls $d/epinano/graphmap/tsv/ | wc -l ) 
qsub -cwd -N epinanom -V -pe smp 4 -b y -o $d/epinano/minimap2/logs2/qout/ -e $d/epinano/minimap2/logs2/qerr/ -t 1-${NUMJOBS_m} /users/enovoa/joramirez/scripts/run_epinano_part2.sge $d/epinano/minimap2
qsub -cwd -N epinanog -V -pe smp 4 -b y -o $d/epinano/graphmap/logs2/qout/ -e $d/epinano/graphmap/logs2/qerr/ -t 1-${NUMJOBS_g} /users/enovoa/joramirez/scripts/run_epinano_part2.sge $d/epinano/graphmap
done


