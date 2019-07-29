#!/bin/bash

echo
echo Parsing information:
array=()
names=()
counter=0
for arg in "$@"; do
	if [[ $arg == *.fastq.gz ]]; then
		gzip -d $arg
		arg=${arg%.gz}; fi
		if [[ $arg == *.fastq ]]; then
			counter=$((counter+1))
			awk '{if(NR%4==1) print $1}' $arg | sed -e "s/^@//" > $OUTPUT/FILE_$counter.txt
			n=$( wc -l $OUTPUT/FILE_$counter.txt | cut -d' ' -f1 )
			sort $OUTPUT/FILE_$counter.txt > $OUTPUT/SORTED_FILE_$counter.txt
			rm $OUTPUT/FILE_$counter.txt
			array=("${array[@]}" "$arg")
			echo input file $counter: $arg contains $n reads
		elif [[ -d $arg ]]; then
			OUTPUT=$arg
			echo The output directory is called $OUTPUT
		else
			names=("${names[@]}" "$arg")
		fi
done


for i in `seq 1 $counter`; do
	echo The input file $i will be named: ${names[i-1]}
	mv $OUTPUT/SORTED_FILE_$i.txt $OUTPUT/SORTED_FILE_${names[i-1]}.txt; done


# Number of Common Reads
echo 
echo Computing number of common reads...

for i in `seq 1 $counter`; do 
	for j in `seq $i $counter`; do 
		if [ $i != $j ]; then
			comm -12 $OUTPUT/SORTED_FILE_${names[i-1]}.txt $OUTPUT/SORTED_FILE_${names[j-1]}.txt > $OUTPUT/COMMON_READS_${names[i-1]}_${names[j-1]}.txt
			nc=$( wc -l $OUTPUT/COMMON_READS_${names[i-1]}_${names[j-1]}.txt | cut -d' ' -f1 )
			echo Common reads in ${names[i-1]} and ${names[j-1]}: $nc
			comm -23 $OUTPUT/SORTED_FILE_${names[i-1]}.txt $OUTPUT/SORTED_FILE_${names[j-1]}.txt > $OUTPUT/READS_in_${names[i-1]}_but_not_in_${names[j-1]}.txt
			nc=$( wc -l $OUTPUT/READS_in_${names[i-1]}_but_not_in_${names[j-1]}.txt | cut -d' ' -f1 )
			echo Reads present in ${names[i-1]} but not in ${names[j-1]}: $nc
			comm -13 $OUTPUT/SORTED_FILE_${names[i-1]}.txt $OUTPUT/SORTED_FILE_${names[j-1]}.txt > $OUTPUT/READS_in_${names[j-1]}_but_not_in_${names[i-1]}.txt
			nc=$( wc -l $OUTPUT/READS_in_${names[j-1]}_but_not_in_${names[i-1]}.txt | cut -d' ' -f1 )
			echo Reads present ${names[j-1]} but not in ${names[i-1]}: $nc
fi; done; done



#per read and quality

echo
echo Computing per-read and quality stats...

for i in `seq 1 $counter`; do 
	echo read_id$'\t'read_length$'\t'mean_quality > $OUTPUT/OUTPUT_${names[i-1]}
	cat $OUTPUT/SORTED_FILE_${names[i-1]}.txt | while read read_id; do
		l=$( grep -A1 $read_id "${array[$i-1]}" | tail -n1 | wc -c )
		grep -A3 $read_id "${array[$i-1]}" | tail -n1 > $OUTPUT/quality
		q=$( python3 scripts/mean_qual.py $OUTPUT/quality )
		echo $read_id$'\t'$l$'\t'$q >> $OUTPUT/OUTPUT_${names[i-1]}
done
echo Computations finished for ${names[i-1]}
done
rm $OUTPUT/quality

echo
echo Plotting

string=""
for i in "${names[@]}"; do string="${string},$i"; done
string="${string:1}"

Rscript scripts/per_read.R -i $OUTPUT -n $string

echo Everything finished without problems

















