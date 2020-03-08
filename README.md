# Best Practices for the Analysis of Oxford Nanopore Direct RNA Sequencing Data

This repository contains scripts to analyze the output of base-called direct RNA sequencing reads, including:
- 1) base-calling quality assessment
- 2) mapping quality assessment
- 3) RNA modification analyses
- 4) PolyA tail estimations

The analysis is thought as a comparison (e.g. comparison of different base-callers, comparison of different mappers, etc)

## Pre-requisites

The following software and modules have been used:

	
|Software|R Modules|
|--|--|
|<table> <tr><th>Software</th><th>Version</th></tr><tr><td>python</td><td>3.6.4</td></tr><tr><td>ont_fast5_api</td><td>1.4.1</td></tr><tr><td>albacore</td><td>2.1.7</td></tr><tr><td>albacore</td><td>2.3.4</td></tr><tr><td>guppy</td><td>2.3.1</td></tr><tr><td>guppy</td><td>3.0.3</td></tr><tr><td>minimap2</td><td>2.16-r922</td></tr><tr><td>graphmap</td><td>0.5.2</td></tr><tr><td>samtools</td><td>1.9</td></tr><tr><td>EpiNano</td><td>1.1</td></tr><tr><td>nanopolish</td><td>0.10.2</td></tr> </table>| <table> <tr><th>Module</th><th>Version</th></tr><tr><td>ggplot2</td><td>3.1.1</td></tr><tr><td>ggExtra</td><td>0.8</td></tr><tr><td>optparse</td><td>1.6.2</td></tr><tr><td>ggpubr</td><td>0.2</td></tr><tr><td>reshape2</td><td>1.4.3</td></tr><tr><td>ggtern</td><td>3.1.0</td></tr><tr><td>lattice</td><td>0.20-38</td></tr><tr><td>latticeExtra</td><td>0.6-28</td></tr><tr><td>KernSmooth</td><td>2.23</td></tr><tr><td>tailfindr</td><td>0.1.0</td></tr>  </table>|

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

* Comparison of sorted.bam files
```
./mapping_analysis_comparison.sh -i ${INPUT_DIR} -o ${OUTPUT_DIR} -r ${REFERENCE_FASTA} -n ${NAMES} ${MORE_OPTIONAL_PARAMETERS} #Where input directory contains the sorted.bam files. It outputs one .csv file per base-caller
#example 1: ./mapping_analysis_comparison.sh -i example_data/ -o output/ -r example_data/reference.fasta -n "graphmap,minimap2"
#example 2: ./mapping_analysis_comparison.sh -i example_data/ -o output/ -r example_data/reference.fasta -n "graphmap,minimap2" -m "A" -t "my_title" -k 0.1 -z 0.8
```
If -m ${base} option is included, extra filtered .csv files will be created containing the mismatch information for the 5-mers that only contain that given base in the central position.
This will automatically return boxplots per base-caller, barplots per bases, boxplots per bases containing mismatch information, and some default ternary plots by base-caller. If further tune of parameters is desired, the newxt command line can be executed:

```
Rscript scripts/mismatch.R -i ${INPUT_DIR} -n ${NAMES} ${MORE_OPTIONAL_PARAMETERS}
#Where ${INPUT_DIR} contains .STATS and .mimatch files (output from previous step) and ${NAMES} would contain the base-caller names separeted by commas
optional arguments:
	-m
		base [A, C, G, T] used as modified for considering 5-mers with m as central position and computing the mismatch pattern and the ternary diagrams.
	-t
		optional plot title
	-k
		threshold for removing positions with lower coverage than this percentage. Default = 0.1
	-z
		zoom for plotting the mismatch pattern. Default = 0.8
	-e
		BOOLEAN for studying the epinano output
#example 1: Rscript mismatch.R -i ~/{example_data} -n "AL_2.1.7,AL_2.3.4,GU_2.3.1,GU_3.0.3" -m A -e 

```

If -z ${zoom} option is present, an extra plot will be created for the mismatch signature with this value as lower boundary for the y axis.

```
#For comparing modifications with ternary plots:
Rscript scripts/ternary.R -m ${mapper} -b ${basecaller} #It still needs improvement

#For checking replicability in two datasets:
Rscript scripts/replicability -i ${INPUT_OUTPUT_DIRECTORY} -n ${NAMES} #The input directory contains the output of ./mapping_analysis_comparison (.mismatches and .STATS output)

```

### Step 5: RNA modification analysis

* EpiNano: https://github.com/enovoa/EpiNano  
We use epinano to get per_site information, it produces as output a per_site.var.csv.slided.onekmer.oneline.5mer.csv file

* Comparison between two .csv files (By now, it only accepts two, either unm vs mod or replicates)
```
./modification_analysis_comparison.sh -b ${MODIFIED_BASE} -u ${UNMODIFIED_CSV_FILE} -m ${MODIFIED_CSV_FILE} -e ${BOOLEAN_FOR_BUILDING_MODEL} -n ${DATASET_NAMES}
#example 1: ./modification_analysis_comparison.sh -b "A" -u example_data/unm_per_site.var.csv.slided.onekmer.oneline.5mer.csv -m example_data/m6A_per_site.var.csv.slided.onekmer.oneline.5mer.csv -o output/ -e false -n "UNM,m6A"
```

### Step 6: PolyA length estimation

Pre-processing
```
samtools view -bF 260 ${prebam} > ${bam} #to keep only mapped reads and remove secondary alignments
bedtools intersect -a ${bed} -b ${bam} -wb | cut -f 9,5 | awk '{print $2 " " $1}' | tr ' ' ',' > intersect.bed
```

* Nanopolish: 
```
nanopolish index -d ${RAW_fast5} ${fastq}  #optional -s for sequencing summary or -f if list of many sequencing summaries
nanopolish polya -t ${threads} -r ${fastq} -b ${bam} -g ${ref} > ${nanopolish_output}
```

* Tailfindr:
```
Rscript scripts/tailfindr_polya.R -i ${base-called_fast5s} -o ${output_dir} -n ${tailfindr_output_name}
```

Post-processing
```
./post_polya.sh ${tailfindr_output_name}.tsv ${nanopolish_output} ${intersect.bed} ${tailfindr_output_gene_names} ${nanopolish_output_gene_names}
```
* Analysis
In the next script we compare both tools using three replicates of S. cerevisiae WT and
three of ime4∆ as input data. We also look for differentially expressed gene (WT vs KO) containg kown m6A modifications (info in sites.bed). We will be updating this script for more general cases.
```
Rscript polyA_guppy.R sites.bed ${tailfindr_DIR} ${nanopolish_DIR}
```

## Citation

We are currently preparing a draft for bioRxiv. 


