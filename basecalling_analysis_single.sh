#This script outputs the number of base-called reads from a given output fastq.gz files 

#FASTQ_FILE=$1

#n=$( awk '{if(NR%4==1) print $1}' ${FASTQ_FILE} | sed -e "s/^@//" | wc -l )
#echo There are $n reads in the output file

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
