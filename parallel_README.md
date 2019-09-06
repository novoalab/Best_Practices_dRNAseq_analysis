### Pipeline for parallel execution:

By now, it is thought to be executed in ant-cluster (CRG cluster)

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


./parallelcripts//pipeline_part4_epinano_clean_part1.sh $INPUT $ref #This command will submit a couple of jobs

./parallelcripts//pipeline_part4_epinano_clean_part2.sh $INPUT #This command will submit a couple of jobs



qlogin -q interactive #We are entering in a cluster node

#Now, we define the variable INPUT [INPUT=/whatever/path]

for d in ${INPUT}/ALBACORE_2.1.7_OUTPUT ${INPUT}/ALBACORE_2.3.4_OUTPUT ${INPUT}/GUPPY_2.3.1_OUTPUT ${INPUT}/GUPPY_3.0.3_OUTPUT; do
cd $d/epinano/minimap2/var.freq
cat *.freq | python2 ~/Epinano/combine_multi_per_site_var_frq.py > $d/epinano/minimap2/per_site.var.csv
python2 ~/Epinano/slide_per_site_var.py $d/epinano/minimap2/per_site.var.csv 5 > $d/epinano/minimap2/per_site.var.sliding.win.csv 
cd $d/epinano/graphmap/var.freq
cat *.freq | python2 ~/Epinano/combine_multi_per_site_var_frq.py > $d/epinano/graphmap/per_site.var.csv
python2 ~/Epinano/slide_per_site_var.py $d/epinano/graphmap/per_site.var.csv 5 > $d/epinano/graphmap/per_site.var.sliding.win.csv 
done

exit



#Now, we have obtained the epinano output, which contains information of all 5-mers. For the curlcakes, we need to filter these kmers for keeping the kmers that only contain the modified base in its center position

./parallelcripts//pipeline_part4_epinano_clean_part3.sh $INPUT $base $name # where base is the base in which the modification takes place (A/C/G/T) and name is the name of the modification (hm5C, m5C...) [base=whateverbase   name=whatevermodification]. If the dataset in INPUT is unmodified (a control) we should also run this script for the base we are interested in [e.g.: base=C name=UNM_C]



#Next step is traingin the SVM, coming soon :)
```
Now, we can run STEPS 2, 4 and 5 from README.md
