# Best Practices for the Analysis of Oxford Nanopore Direct RNA Sequencing Data

This repository contains scripts to analyze the output of base-called direct RNA sequencing reads, including:
- 1) base-calling quality assessment
- 2) mapping quality assessment
- 3) RNA modification analyses
- 4) PolyA tail estimations

The analysis can be ran for a single fastq file, of in a comparative mode (e.g. comparison of different base-callers, comparison of different mappers, etc)

## Pre-requisites

The following software and modules have been used:

	
|Software|R Modules|
|--|--|
|<table> <tr><th>Software</th><th>Version</th></tr><tr><td>python</td><td>3.6.4</td></tr><tr><td>ont_fast5_api</td><td>1.4.1</td></tr><tr><td>albacore</td><td>2.1.7</td></tr><tr><td>albacore</td><td>2.3.4</td></tr><tr><td>guppy</td><td>2.3.1</td></tr><tr><td>guppy</td><td>3.0.3</td></tr><tr><td>biopython</td><td>1.73</td></tr><tr><td>minimap2</td><td>2.16-r922</td></tr><tr><td>graphmap</td><td>0.5.2</td></tr><tr><td>samtools</td><td>1.9</td></tr><tr><td>EpiNano</td><td>1.1</td></tr> </table>| <table> <tr><th>Module</th><th>Version</th></tr><tr><td>ggplot2</td><td>3.1.1</td></tr><tr><td>ggExtra</td><td>0.8</td></tr><tr><td>optparse</td><td>1.6.2</td></tr><tr><td>ggpubr</td><td>0.2</td></tr><tr><td>reshape2</td><td>1.4.3</td></tr><tr><td>ggtern</td><td>3.1.0</td></tr><tr><td>lattice</td><td>0.20-38</td></tr><tr><td>latticeExtra</td><td>0.6-28</td></tr> </table>|

## Analysis Steps 

### Step 1: Base-calling

Albacore doesn't support multifast5 files, if these are the input files, conversion from single_fast5 to multi_fast5 is needed:
>
	multi_to_single_fast5 -i ${folder_with_multifast5_files} -s ${output_folder} -t ${threads_number} 
	
* Albacore v2.1.7 & v2.3.4:
```{r, engine='bash', count_lines}
read_fast5_basecaller.py --flowcell ${FLOWCELL} --kit ${KIT} --output_format fastq,fast5 -n ${NUMFAST5} --input ${INPUT_DIRECTORY} --save_path ${OUTPUT_DIRECTORY} --worker_threads ${NUMBER_OF_THREADS} --disable_filtering
```

* Guppy v2.3.1 & v3.0.3:
```{r, engine='bash', count_lines}
guppy_basecaller --flowcell ${FLOWCELL} --kit ${KIT} --fast5_out --input ${INPUT_DIRECTORY} --save_path ${OUTPUT_DIRECTORY} --cpu_threads_per_caller ${NUMBER_OF_THREADS}
```
### Step 2: Analysis of base-calling

* Single fastq file
```{r}
./basecalling_analysis_single.sh ${FASTQ_FILE}
#example: ./basecalling_analysis_single.sh example_data/test_1.fastq
```
* Comparison of fastq files (the second part of this script is thought to compare base-callers as read_ids are used)
```{r}
./basecalling_analysis_comparison.sh ${OUTPUT_DIRECTORY} ${ALL_FASTQ_FILES} ${ALL_FASTQ_NAMES_FOR_PLOTTING}
#example: ./basecalling_analysis_comparison.sh output/ example_data/test_1.fastq example_data/test_1.fastq dataset_1 dataset_2
```

### Step 3: Mapping

'U' to 'T' conversion:
```{r, engine='bash', count_lines}
for i in *fastq; do
awk '{ if (NR%4 == 2) {gsub(/U/,"T",$1); print $1} else print }' $i > ${i%.fastq}.U2T.fastq;
rm $i
done
```
* minimap2 default:
```{r, engine='bash', count_lines}
minimap2 -ax map-ont ${FASTA_REFERENCE} ${FASTQ_FILE} > ${FASTQ_FILE%.U2T*}.sam
```
* minimap2 sensitive:
```
minimap2 -ax map-ont -k 5 ${FASTA_REFERENCE} ${FASTQ_FILE} > ${FASTQ_FILE%.U2T*}.sam   #-w5 and -m20 are also good for increasing sensitivity
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


### Step 4: Analysis of mapping

* Single sorted.bam file
```
```

* Comparison of sorted.bam files
```
./mapping_analysis_comparison.sh ${DIRECTORY_EITH_SORTED.BAM_FILES} ${OUTPUT_DIR} ${REFERENCE_FASTA} ${NAMES} ${MORE_OPTIONAL_PARAMETERS}
#example 1: ./mapping_analysis_comparison.sh example_data/ output/ example_data/reference.fasta "graphmap,minimap2"
#example 2: ./mapping_analysis_comparison.sh example_data/ output/ example_data/reference.fasta "graphmap,minimap2" "A" "my_title" 0.1 0.8

#For comparing modifications with ternary plots:
Rscript ternary.R -m ${mapper} -b ${basecaller} #It still needs improvement
```

### Step 5: RNA modification analysis

* EpiNano: https://github.com/enovoa/EpiNano  
We use epinano to get per_site information, it produces as output a per_site.var.csv.slided.onekmer.oneline.5mer.csv file and we then have to filter these results:
```
base=T  #base is going to be the nucleotide we want to study depending on the RNA modification
head -1 per_site.var.csv.slided.onekmer.oneline.5mer.csv > per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv
awk -F"," "\$1 ~ /[^/$base/][^/$base/][/$base/][^/$base/][^/$base/]/" per_site.var.csv.slided.onekmer.oneline.5mer.csv >> per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv
mv per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv ${name}_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv #we change here the name according to the specific modification
```
We then divide our data into training and test:
```
cut -d "," -f 5-14,20-24 ${unm}_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv > test
tail -n +2 test > test2
sed 's/$/,UNM/' test2 > test3
a=$( wc -l test3 )
b=$( echo $a | awk '{print $1;}' )
t=$( echo "$b * 0.75" | bc -l | xargs printf "%.*f\n" 0 )   # 75 % for training
p=$( echo "$b * 0.25" | bc -l | xargs printf "%.*f\n" 0 )   # 25 % for predicting
echo q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,del1,del2,del3,del4,del5,sample > training
echo q1,q2,q3,q4,q5,mis1,mis2,mis3,mis4,mis5,del1,del2,del3,del4,del5,sample > predicting
head -$t test3 >> training 
tail -$p test3 >> predicting

cut -d "," -f 5-14,20-24 ${mod}_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv > test
tail -n +2 test > test2
sed 's/$/,MOD/' test2 > test3
head -$t test3 >> training 
tail -$p test3 >> predicting
rm test test2 test3 
```
And we finally run Epinano's SVM script:
```
SVM.py -a -t training -p predicting -cl 1-15 -mc 16 -o $mod
```

* Analysis
```
Rscript kmer.R -i ${input_directory} #The input directory would contain ${unm}_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv and ${mod}_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv
Rscript PCA.R -b ${basecaller} -m ${mapper} #This script takes as input the name of the base-caller and mapper used, names that will be used for the output directory containing the PCA plots. 
```

## Citation

We are currently preparing a draft for bioRxiv. 


