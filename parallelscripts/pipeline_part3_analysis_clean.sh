#!/bin/bash

INPUT=$1
ref=$2

###Number of sequenced reads

n=$( find $INPUT/singleRAW/ -type f | wc -l )
num=$(($n - 1))
echo $num reads were sequenced


###Number of base-called reads

for d in ${INPUT}/ALBACORE_2.1.7_OUTPUT ${INPUT}/ALBACORE_2.3.4_OUTPUT ${INPUT}/GUPPY_2.3.1_OUTPUT ${INPUT}/GUPPY_3.0.3_OUTPUT; do
fastq_file=${d}/fastq/RNA*
gzip -d $fastq_file
awk '{if(NR%4==1) print $1}' ${fastq_file%.gz} | sed -e "s/^@//" > ${d}/fastq/basecalled_nonsorted_reads
sort ${d}/fastq/basecalled_nonsorted_reads > ${d}/fastq/basecalled_reads
rm ${d}/fastq/basecalled_nonsorted_reads
n=$( cat ${d}/fastq/basecalled_reads | wc -l )
s=${d%_*}
echo $n reads were basecalled using ${s,,}
done


### Common reads

#Changing location and name
for d in ALBACORE_2.1.7_OUTPUT ALBACORE_2.3.4_OUTPUT GUPPY_2.3.1_OUTPUT GUPPY_3.0.3_OUTPUT; do
mkdir -p $INPUT/COMPARISONS
file=${INPUT}/${d}/fastq/basecalled_reads
if [[ $d == ALBACORE_2.1.7_OUTPUT/ ]]
then
mv $file $INPUT/COMPARISONS/ALBACORE_O
elif [[ $d == ALBACORE_2.3.4_OUTPUT/ ]]
then
mv $file $INPUT/COMPARISONS/ALBACORE_N
elif [[ $d == GUPPY_2.3.1_OUTPUT/ ]]
then
mv $file $INPUT/COMPARISONS/GUPPY_O
elif [[ $d == GUPPY_3.0.3_OUTPUT/ ]]
then
mv $file $INPUT/COMPARISONS/GUPPY_N
else
echo wrong input
fi
done


C=$INPUT/COMPARISONS
comm -12 $C/ALBACORE_N $C/ALBACORE_O > $C/ALBACORE_COMMON
comm -13 $C/ALBACORE_N $C/ALBACORE_O > $C/ALBACORE_O_ONLY
comm -23 $C/ALBACORE_N $C/ALBACORE_O > $C/ALBACORE_N_ONLY
comm -12 $C/GUPPY_N $C/GUPPY_O > $C/GUPPY_COMMON
comm -13 $C/GUPPY_N $C/GUPPY_O > $C/GUPPY_O_ONLY
comm -23 $C/GUPPY_N $C/GUPPY_O > $C/GUPPY_N_ONLY


cd $C
touch basecalled_information.txt
a=$( cat $C/ALBACORE_COMMON | wc -l )
echo Both versions of Albacore base-called $a reads in common >> basecalled_information.txt
echo Both versions of Albacore base-called $a reads in common
a=$( cat $C/ALBACORE_O_ONLY | wc -l )
echo Albacore 2.1.7 uniquely base-called $a reads >> basecalled_information.txt
echo Albacore 2.1.7 uniquely base-called $a reads
a=$( cat $C/ALBACORE_N_ONLY | wc -l )
echo Albacore 2.3.4 uniquely base-called $a reads >> basecalled_information.txt
echo Albacore 2.3.4 uniquely base-called $a reads

a=$( cat $C/GUPPY_COMMON | wc -l )
echo Both versions of Guppy base-called $a reads in common >> basecalled_information.txt
echo Both versions of Guppy base-called $a reads in common
a=$( cat $C/GUPPY_O_ONLY | wc -l )
echo Guppy 2.3.4 uniquely base-called $a reads >> basecalled_information.txt
echo Guppy 2.3.4 uniquely base-called $a reads
a=$( cat $C/GUPPY_N_ONLY | wc -l )
echo Guppy 3.0.3 uniquely base-called $a reads >> basecalled_information.txt
echo Guppy 3.0.3 uniquely base-called $a reads

for i in $C/ALBACORE_COMMON $C/ALBACORE_O_ONLY $C/ALBACORE_N_ONLY $C/GUPPY_COMMON $C/GUPPY_O_ONLY $C/GUPPY_N_ONLY; do
rm $i
done


###Per read


mkdir $C/INPUT

cd $C
for i in ALBACORE_O ALBACORE_N GUPPY_O GUPPY_N
do
	split -l 500 $i INPUT/${i}_
	rm $i
done

mkdir -p $C/OUTPUT $C/logs/qout $C/logs/qerr
NUMFILES=$( find $C/INPUT/ -type f | wc -l )

AL_O=$INPUT/ALBACORE_2.1.7_OUTPUT/fastq/RNA*
AL_N=$INPUT/ALBACORE_2.3.4_OUTPUT/fastq/RNA*
GU_O=$INPUT/GUPPY_2.3.1_OUTPUT/fastq/RNA*
GU_N=$INPUT/GUPPY_3.0.3_OUTPUT/fastq/RNA*

qsub -cwd -N read -V -t 1-${NUMFILES} -o $C/logs/qout/ -e $C/logs/qerr/ /users/enovoa/joramirez/scripts/computation_parallel.sh $AL_O $AL_N $GU_O $GU_N $C/ $C/OUTPUT


rm $C/quality_*
for i in $AL_O $AL_N $GU_O $GU_N; do
gzip $i
done

echo read_id$'\t'read_length$'\t'mean_quality > $C/OUTPUT/ALBACORE_N_OUTPUT
for i in $C/OUTPUT/ALBACORE_N_*_OUTPUT  
do cat $i >> $C/OUTPUT/ALBACORE_N_OUTPUT
rm $i
done
echo read_id$'\t'read_length$'\t'mean_quality > $C/OUTPUT/ALBACORE_O_OUTPUT
for i in $C/OUTPUT/ALBACORE_O_*_OUTPUT  
do cat $i >> $C/OUTPUT/ALBACORE_O_OUTPUT
rm $i
done
echo read_id$'\t'read_length$'\t'mean_quality > $C/OUTPUT/GUPPY_N_OUTPUT
for i in $C/OUTPUT/GUPPY_N_*_OUTPUT    
do cat $i >> $C/OUTPUT/GUPPY_N_OUTPUT
rm $i
done
echo read_id$'\t'read_length$'\t'mean_quality > $C/OUTPUT/GUPPY_O_OUTPUT
for i in $C/OUTPUT/GUPPY_O_*_OUTPUT    
do cat $i >> $C/OUTPUT/GUPPY_O_OUTPUT
rm $i
done


###Number of mapped reads

for d in ${INPUT}/ALBACORE_2.1.7_OUTPUT ${INPUT}/ALBACORE_2.3.4_OUTPUT ${INPUT}/GUPPY_2.3.1_OUTPUT ${INPUT}/GUPPY_3.0.3_OUTPUT; do
for i in $d/bam/*.sorted.bam; do
echo $i
samtools flagstat $i; done; done


###From bam to stats:


for d in ${INPUT}/ALBACORE_2.1.7_OUTPUT ${INPUT}/ALBACORE_2.3.4_OUTPUT ${INPUT}/GUPPY_2.3.1_OUTPUT ${INPUT}/GUPPY_3.0.3_OUTPUT; do
cd $d/bam/
for i in *.sorted.bam; do
samtools mpileup -f $ref $i -o $i.mpileup   #remember that now the reference for 71 % es diferente.
done
/users/enovoa/joramirez/scripts/from_bam_to_stats.sh
done





