for i in *sorted.bam; do
	# Convert qualities to ordinal numbers
	python scripts/1.base_quality_chr2ord.keep_all.py $i.mpileup > $i.mpileup.digit

	# Get mean and median quality scores (produces output with header)
	python scripts/get_mean_median_qualities.py $i.mpileup.digit > $i.mean_median_quals

	# Get general file info (coverage, ref_nuc)
	echo -ne "chr\tpos\tref_nuc\tcoverage\n" > $i.info
	cut -f 1-4 $i.mpileup >> $i.info

	# Get number of rt-stop, insertion and deletions (produces file with header)
	echo "0" > $i.mpileup.digit.rtstop.tmp
	cat $i.mpileup.digit | awk {'print $5'} | awk -F "$" '{print NF-1}' >> $i.mpileup.digit.rtstop.tmp 	 #NF = number of fields in the record
	cat $i.mpileup.digit | awk {'print $5'} | awk -F "-" '{print NF-1}'  > $i.mpileup.digit.del
	cat $i.mpileup.digit | awk {'print $5'} | awk -F "+" '{print NF-1}'  > $i.mpileup.digit.ins
	sed '$ d' $i.mpileup.digit.rtstop.tmp > $i.mpileup.digit.rtstop # remove last line as we pushed everything 1
	echo -ne "rtstop\tins\tdel\n" > $i.3moreinfo
	paste $i.mpileup.digit.rtstop $i.mpileup.digit.del $i.mpileup.digit.ins >> $i.3moreinfo
	
	# Get mismatches
	perl scripts/pileup2base-master/pileup2base.pl  $i.mpileup 0 $i.mismatches.tmp 
	cut -f 4-12 $i.mismatches.tmp | tail -n +2 > $i.mismatches.tmp2
	cut -f 1-3 $i.mismatches.tmp | tail -n +2 > $i.mismatches.tmp3
	echo -ne "chr\tpos\tref_nuc\tA\tT\tC\tG\n" > $i.mismatches
 	awk 'BEGIN{print "sum" > "kk"}NR>0{print $1+$5"\t"$2+$6"\t"$3+$7"\t"$4+$8}' $i.mismatches.tmp2 > $i.mismatches.tmp4
	paste $i.mismatches.tmp3 $i.mismatches.tmp4 >> $i.mismatches
	rm $i.mismatches.tmp $i.mismatches.tmp2 $i.mismatches.tmp3 $i.mismatches.tmp4 kk
 
	# Identity of insertions and deletions
	perl scripts/pileup2base-master/pileup2baseindel.pl -i $i.mpileup
	cut -f 12-13 sample1.txt > $i.indel_info

	# Merge mismatches and indel info
	paste $i.mismatches $i.indel_info > $i.mismatch_and_indels

	# Merge all files
	paste $i.info $i.mean_median_quals $i.3moreinfo > $i.STATS

	# Remove temporary files
	rm $i.info $i.mean_median_quals $i.3moreinfo $i.mpileup.digit.rtstop $i.mpileup.digit.del $i.mpileup.digit.ins $i.mpileup.digit.rtstop.tmp $i.mpileup.digit sample1.txt

done
