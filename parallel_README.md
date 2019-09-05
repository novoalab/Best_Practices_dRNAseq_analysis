*** Pipeline for parallel execution:

```
mkdir ${OUTPUT}
mkdir ${OUTPUT}/singleRAW

qlogin
multi_to_single_fast5 -i $RAW -s ${OUTPUT}/singleRAW/ -t ${n_threads}   # RAW is a path of the directory with raw fast5s
exit

./parallelscripts/pipeline_part1_basecalling.sh $RAW $OUTPUT $multifast5   # RAW is a directory with the raw data; OUTPUT, the name of the directory that will contain the output files, and multifast5, a boolean for specifying if the input raw data is multifast5 or single fast5 [if the data is multifast5, multifast5=true]
#e.g.: ./pipeline_part1_basecalling.sh ~/RAW/RNA936935_fragmented/fast5/ ~/RNA936935_fragmented true 

./parallelcripts/pipeline_part2_mapping.sh $INPUT $ref   # INPUT is the $OUTPUT of the previous step; ref is the reference that will be used to map

./parallelcripts/pipeline_part3_analysis_clean.sh $INPUT $ref   # INPUT is the $OUTPUT of the previous step; ref is the reference that will be used to map
```
