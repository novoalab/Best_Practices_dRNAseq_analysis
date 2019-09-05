#!/bin/bash
INPUT=$1      # directory containing four folders: one per base-caller
ref=$2


###Combining fastq files

for d in ALBACORE_2.1.7_OUTPUT ALBACORE_2.3.4_OUTPUT GUPPY_2.3.1_OUTPUT GUPPY_3.0.3_OUTPUT; do
b=${INPUT%\/}_${a,,}; 
x=${b##*/};
cat ${INPUT}/${d}/fastq/fastq_runid* > ${INPUT}/${d}/fastq/${x}.fastq.gz;
#rm $d/fastq/fastq_runid*;
done

for d in ALBACORE_2.1.7_OUTPUT ALBACORE_2.3.4_OUTPUT GUPPY_2.3.1_OUTPUT GUPPY_3.0.3_OUTPUT; do b=${INPUT%\/}_${a,,}; x=${b##*/}; echo $x; cat ${INPUT}/${d}/fastq/fastq_runid* > ${b}.fastq.gz; echo ${INPUT}/${d}/${x}; done


###Mapping

#From U to T

for d in ${INPUT}/ALBACORE_2.1.7_OUTPUT ${INPUT}/ALBACORE_2.3.4_OUTPUT ${INPUT}/GUPPY_2.3.1_OUTPUT ${INPUT}/GUPPY_3.0.3_OUTPUT; do
fastq_file=${d}/fastq/RNA*
name=$( basename ${fastq_file} )
gzip -d $fastq_file
awk '{ if (NR%4 == 2) {gsub(/U/,"T",$1); print $1} else print }' ${fastq_file%/*}/${name%.gz} > ${fastq_file%/*}/${name%.fastq}_UtoT.fastq
gzip ${fastq_file%/*}/${name%.fastq}_UtoT.fastq
rm ${fastq_file%/*}/${name%.gz}
done


for d in ${INPUT}/ALBACORE_2.1.7_OUTPUT ${INPUT}/ALBACORE_2.3.4_OUTPUT ${INPUT}/GUPPY_2.3.1_OUTPUT ${INPUT}/GUPPY_3.0.3_OUTPUT; do
fastq_file=${d}/fastq/RNA*
sam_name=$( basename ${fastq_file} )
mkdir -p $d/bam
qsub -cwd -N minimap2d -V -o $d/logs/qout -e $d/logs/qerr /users/enovoa/joramirez/scripts/mapping.sge $ref $fastq_file $d/bam/${sam_name%fastq.gz}_minimap2_default.sam minimap2_default 
qsub -cwd -N minimap2s -V -o $d/logs/qout -e $d/logs/qerr /users/enovoa/joramirez/scripts/mapping.sge $ref $fastq_file $d/bam/${sam_name%fastq.gz}_minimap2_sensitive.sam minimap2_sensitive 
qsub -cwd -N graphmapd -V -o $d/logs/qout -e $d/logs/qerr /users/enovoa/joramirez/scripts/mapping.sge $ref $fastq_file $d/bam/${sam_name%fastq.gz}_graphmap_default.sam graphmap_default 
qsub -cwd -N graphmaps -V -o $d/logs/qout -e $d/logs/qerr /users/enovoa/joramirez/scripts/mapping.sge $ref $fastq_file $d/bam/${sam_name%fastq.gz}_graphmap_sensitive.sam graphmap_sensitive 
done




