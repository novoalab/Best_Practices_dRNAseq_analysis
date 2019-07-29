#!/bin/bash
INPUT_DIR=$1 #folder with *.sorted.bam files
OUTPUT_DIR=$2
ref=$3 #reference
names=$4 #string with the dataset names separated with commas: e.g.: "1,2"

if [ -z "$5" ]; then
	mod_log=false; else
	mod_log=true
	mod=$5;fi

if [ -z "$6" ]; then
	title_log=false; else
	title_log=true
	title=$6;fi

if [ -z "$7" ]; then
	threshold=0.1; else
	threshold=$7;fi

if [ -z "$8" ]; then
	zoom=0.8; else
	zoom=$8;fi

for i in $INPUT_DIR/*sorted.bam; do
	# Compute the mpileup
	samtools mpileup -f $ref $i -o $i.mpileup
	OUT=$OUTPUT_DIR/${i##*/}
	# Convert qualities to ordinal numbers
	python scripts/1.base_quality_chr2ord.keep_all.py $i.mpileup > $OUT.mpileup.digit

	# Get mean and median quality scores (produces output with header)
	python scripts/get_mean_median_qualities.py $OUT.mpileup.digit > $OUT.mean_median_quals

	# Get general file info (coverage, ref_nuc)
	echo -ne "chr\tpos\tref_nuc\tcoverage\n" > $OUT.info
	cut -f 1-4 $i.mpileup >> $OUT.info

	# Get number of rt-stop, insertion and deletions (produces file with header)
	echo "0" > $OUT.mpileup.digit.rtstop.tmp
	cat $OUT.mpileup.digit | awk {'print $5'} | awk -F "$" '{print NF-1}' >> $OUT.mpileup.digit.rtstop.tmp 	 #NF = number of fields in the record
	cat $OUT.mpileup.digit | awk {'print $5'} | awk -F "-" '{print NF-1}'  > $OUT.mpileup.digit.del
	cat $OUT.mpileup.digit | awk {'print $5'} | awk -F "+" '{print NF-1}'  > $OUT.mpileup.digit.ins
	sed '$ d' $OUT.mpileup.digit.rtstop.tmp > $OUT.mpileup.digit.rtstop # remove last line as we pushed everything 1
	echo -ne "rtstop\tins\tdel\n" > $OUT.3moreinfo
	paste $OUT.mpileup.digit.rtstop $OUT.mpileup.digit.del $OUT.mpileup.digit.ins >> $OUT.3moreinfo
	
	# Get mismatches
	perl scripts/pileup2base.pl  $i.mpileup 0 $OUT.mismatches.tmp 
	cut -f 4-12 $OUT.mismatches.tmp | tail -n +2 > $OUT.mismatches.tmp2
	cut -f 1-3 $OUT.mismatches.tmp | tail -n +2 > $OUT.mismatches.tmp3
	echo -ne "chr\tpos\tref_nuc\tA\tT\tC\tG\n" > $OUT.mismatches
 	awk 'BEGIN{print "sum" > "kk"}NR>0{print $1+$5"\t"$2+$6"\t"$3+$7"\t"$4+$8}' $OUT.mismatches.tmp2 > $OUT.mismatches.tmp4
	paste $OUT.mismatches.tmp3 $OUT.mismatches.tmp4 >> $OUT.mismatches
	rm $OUT.mismatches.tmp $OUT.mismatches.tmp2 $OUT.mismatches.tmp3 $OUT.mismatches.tmp4 kk

	# Merge all files
	paste $OUT.info $OUT.mean_median_quals $OUT.3moreinfo > $OUT.STATS

	# Remove temporary files
	rm $OUT.info $OUT.mean_median_quals $OUT.3moreinfo $OUT.mpileup.digit.rtstop $OUT.mpileup.digit.del $OUT.mpileup.digit.ins $OUT.mpileup.digit.rtstop.tmp $OUT.mpileup.digit

done

if [ $title_log = true ] && [ $mod_log = true ]; then
Rscript scripts/mismatch.R -i $OUTPUT_DIR -n $names -e -m $mod -t $title -k $threshold -z $zoom
elif [ $title_log = true ]; then
Rscript scripts/mismatch.R -i $OUTPUT_DIR -n $names -e -t $title -k $threshold -z $zoom
elif [ $mod_log = true ]; then
Rscript scripts/mismatch.R -i $OUTPUT_DIR -n $names -e -m $mod -k $threshold -z $zoom;
else
Rscript scripts/mismatch.R -i $OUTPUT_DIR -n $names -e -k $threshold -z $zoom;
fi


#maybe accept sam and do the post-processing here, and also accept bam and do the sorting and indexing here.
#Use names for the first part too maybe


