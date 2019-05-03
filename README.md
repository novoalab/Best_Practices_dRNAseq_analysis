# Best Practices for Direct RNA Sequencing Analysis

### Pre-requisits

|**Software**| **Version** |
|:---------:|-------------|
| python    | 3.6.4 |
| albacore  | 1.2.7    |
| guppy v2 | 2.3.1  |
| guppy v3  | 3.0.3  |
| POREquality | 0.9.8 |
| minimap2  | 2.16-r922    |
| graphmap  | 0.5.2  |


### Step 1: Base-calling

Albacore v2.1.7:
read_fast5_basecaller.py --flowcell ${FLOWCELL} --kit ${KIT} --output_format fastq,fast5 -n ${NUMFAST5} --input ${INPUT_DIRECTORY} --save_path ${OUTPUT_DIRECTORY} --worker_threads ${NUMBER_OF_THREADS} --disable_filtering

Guppy v2.3.1 & 3.0.3:
```{r, engine='bash', count_lines}
guppy_basecaller --flowcell ${FLOWCELL} --kit ${KIT} --fast5_out --input ${INPUT_DIRECTORY} --save_path ${OUTPUT_DIRECTORY} --cpu_threads_per_caller ${NUMBER_OF_THREADS}
```
### Step 2: Analysis of base-calling



