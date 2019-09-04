#!/bin/bash



tail -n +2 $1 > noheader_tailfindr
echo read_id$'\t'gene_name$'\t'tail_start$'\t'tail_end$'\t'samples_per_nt$'\t'tail_length > $4
awk -v OFS='\t' -F',' 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){next}; print $1,a[$1],$2,$3,$4,$5}' $3 noheader.tailfindr >> $4

grep PASS $2 > noheader_nanopolish
echo readname$'\t'contig$'\t'position$'\t'leader_start$'\t'adapter_start$'\t'polya_start$'\t'transcript_start$'\t'read_rate$'\t'polya_length$'\t'gene_name > $5
awk -v OFS='\t' -F'\t' 'FNR==NR{a[$1]=$2;next}{if(a[$1]==""){next}; print $1,$2,$3,$4,$5,$6,$7,$8,$9,a[$1]}' $3 noheader_nanopolish >> $5

rm noheader_tailfindr noheader_nanopolish
