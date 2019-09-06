#!/bin/bash

INPUT=$1
ref=$2

source /users/enovoa/joramirez/.bashrc


for d in ${INPUT}/ALBACORE_2.1.7_OUTPUT ${INPUT}/ALBACORE_2.3.4_OUTPUT ${INPUT}/GUPPY_2.3.1_OUTPUT ${INPUT}/GUPPY_3.0.3_OUTPUT; do
if [[ ! -f $d/epinano/minimap2/logs/qout ]]
then
 mkdir -p $d/epinano/minimap2/logs/qout
 mkdir -p $d/epinano/minimap2/logs/qerr
 mkdir -p $d/epinano/minimap2/bam
 mkdir -p $d/epinano/minimap2/parts
 mkdir -p $d/epinano/minimap2/tsv
 mkdir -p $d/epinano/minimap2/var.freq
 mkdir -p $d/epinano/graphmap/logs/qout
 mkdir -p $d/epinano/graphmap/logs/qerr
 mkdir -p $d/epinano/graphmap/bam
 mkdir -p $d/epinano/graphmap/parts
 mkdir -p $d/epinano/graphmap/tsv
 mkdir -p $d/epinano/graphmap/var.freq
else
 echo "[ WARNING ] Found existing files. Clean up?"
fi
fastq_file=${d}/fastq/RNA*
qsub -cwd -N epinanom -V -pe smp 4 -b y -o $d/epinano/minimap2/logs/qout/ -e $d/epinano/minimap2/logs/qerr/ /users/enovoa/joramirez/scripts/run_epinano.sge $fastq_file $d/epinano/minimap2 $ref minimap2
qsub -cwd -N epinanog -V -pe smp 4 -b y -o $d/epinano/graphmap/logs/qout/ -e $d/epinano/graphmap/logs/qerr/ /users/enovoa/joramirez/scripts/run_epinano.sge $fastq_file $d/epinano/graphmap $ref graphmap
done

