# Best Practices for Direct RNA Sequencing Analysis

### Pre-requisits

The following software and modules have been used:

|**Software**| **Version** |
|:---------:|-------------|
| python    | 3.6.4 |
| albacore v1 | 2.1.7 |
| albacore v2 | 2.3.4 |
| guppy v2 | 2.3.1  |
| guppy v3  | 3.0.3  |
| POREquality | 0.9.8 |
| EMBOSS stretcher | 6.6.0.0 |
| biopython | 1.73 |
| minimap2  | 2.16-r922    |
| graphmap  | 0.5.2  |
| samtools | 1.9 |


### Step 1: Base-calling

* Albacore v2.1.7 & v2.3.4:
```{r, engine='bash', count_lines}
read_fast5_basecaller.py --flowcell ${FLOWCELL} --kit ${KIT} --output_format fastq,fast5 -n ${NUMFAST5} --input ${INPUT_DIRECTORY} --save_path ${OUTPUT_DIRECTORY} --worker_threads ${NUMBER_OF_THREADS} --disable_filtering
```

* Guppy v2.3.1 & v3.0.3:
```{r, engine='bash', count_lines}
guppy_basecaller --flowcell ${FLOWCELL} --kit ${KIT} --fast5_out --input ${INPUT_DIRECTORY} --save_path ${OUTPUT_DIRECTORY} --cpu_threads_per_caller ${NUMBER_OF_THREADS}
```
### Step 2: Analysis of base-calling

* Number of common base-called reads between approaches:
```{r, engine='bash', count_lines}
awk '{if(NR%4==1) print $1}' ${FASTQ_FILE} | sed -e "s/^@//" > ${OUTPUT_FILE}    #It lists all base-called reads from a fastq file
sort ${OUTPUT_FILE} > ${SORTED_OUTPUT_FILE}
wc -l ${SORTED_OUTPUT_FILE}
comm -12 ${SORTED_OUTPUT_FILE_X} ${SORTED_OUTPUT_FILE_Y} > ${COMMON_READS}
comm -13 ${SORTED_OUTPUT_FILE_X} ${SORTED_OUTPUT_FILE_Y} > ${ONLY_PRESENT_IN_Y_READS}
comm -23 ${SORTED_OUTPUT_FILE_X} ${SORTED_OUTPUT_FILE_Y} > ${ONLY_PRESENT_IN_X_READS}
```
* POREquality:
https://github.com/carsweshau/POREquality
It takes as input the sequencing_summary.txt (output from any base-caller), if the base-calling is run in parallel, the sequencing_summary may be created in peaces, if this case we can merge them by:
```{r, engine='bash', count_lines}
header_file=$(ls *txt.gz | shuf -n 1)
zcat $header_file | head -1 > ${FINAL_SEQUENCING_SUMMARY}
for i in sequencing_summary*; do
zcat $i | tail -n +2 >> ${FINAL_SEQUENCING_SUMMARY}
done
```
* Per_read comparison:
```{r, engine='bash', count_lines}
echo $'read_id\talbacore_length\tguppy_length\talignment_length\tmatches\tmismatches\tgaps_albacore\tgaps_guppy\tseq_id\tgaps_albacore_50\tgaps_guppy_50\tinsertions\tdeletions\tinsertions_50\tdeletions_50' > final_output

cat ${COMMON_READS} | while read read_id
do	
	echo '>albacore' > read_A
	grep -A1 $read_id ${FASTQ_FILE_BASECOLLER_A} | tail -n1 >> read_A
	echo '>guppy' > read_B
	grep -A1 $read_id ${FASTQ_FILE_BASECOLLER_B} | tail -n1 >> read_B

	al=$(sed '2q;d' read_A | wc -c)
	gl=$(sed '2q;d' read_B | wc -c)

	stretcher read_A read_B output -aformat_outfile fasta

	a=$(python3 compute_identity.py output)
	echo $a
	echo $read_id$'\t'$al$'\t'$gl$'\t'$a >> final_output
done
rm read_A
rm read_B
```
If it is executed in parallel and we get different final_output_${SGE_TASK_ID}, we can simply merge them:
```{r, engine='bash', count_lines}
for i in final_output_*
do cat $i >> final_output
rm $i
done
```
For plotting the comparison:
```{r, engine='bash', count_lines}
Rscript -e "rmarkdown::render(input='per_read_analysis.Rmd')" -t ${THRESHOLD} -i ${INPUT_FILES} -n ${INPUT_NAMES}
# where ${THRESHOLD} is the upper boundary for removing outliers, ${INPUT_FILES} is a string containing all the input files we want to compare separated by commas, e.g. final_output_RNAAB063141,final_output_RNAAB064293 and ${INPUT_NAMES} is a string containing the names that we want for all the input files, e.g. "0 % DMS,0.5 % DMS"
```

### Step 3: Mapping

'U' to 'T' conversion:
```{r, engine='bash', count_lines}
for i in *fastq; do
awk '{ if (NR%4 == 2) {gsub(/U/,"T",$1); print $1} else print }' $i > ${i%.fastq}.U2T.fastq;
rm $i
done
```
* minimap2:
```{r, engine='bash', count_lines}
minimap2 -ax map-ont ${FASTA_REFERENCE} ${FASTQ_FILE} > ${FASTQ_FILE%.U2T*}.sam
```

* graphmap default:
```{r, engine='bash', count_lines}
graphmap align -r ${FASTA_REFERENCE} -d ${FASTQ_FILE} -o ${FASTQ_FILE%.U2T*}.sam -v 1 -K fastq
```

* graphmap sensitive:
```{r, engine='bash', count_lines}
graphmap align -r ${FASTA_REFERENCE} -d ${FASTQ_FILE} -o ${FASTQ_FILE%.U2T*}.sam -v 1 -K fastq --rebuild-index --double-index --mapq -1 -x sensitive -z -1 --min-read-len 0 -A 7 -k 5
```

Post-processing:
```{r, engine='bash', count_lines}
for i in *.sam; do
str=$(echo $i| sed  's/.sam//')
samtools view -Sbh $i > $str.bam  #Transforming the .sam into .bam
samtools sort -@ 4 $str.bam > $str.sorted.bam #Sorting the bam file
samtools index $str.sorted.bam #Indexing the sorted bam
rm $i #Removing the sam file
rm $str.bam #Removing the not sorted bam
done
```


### Step 4: RNA Modification Analysis

* EpiNano: https://github.com/enovoa/EpiNano  
It produces as output a per_site.var.csv.slided.onekmer.oneline.5mer.csv file and we then filter these results:
```{r, engine='bash', count_lines}
head -1 per_site.var.csv.slided.onekmer.oneline.5mer.csv > per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv
awk -F"," '$1 ~ /[^T][^T][T][^T][^T]/' per_site.var.csv.slided.onekmer.oneline.5mer.csv >> per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv
```


