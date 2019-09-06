#!/bin/bash


INPUT=$1
base=$2
name=$3


for d in ${INPUT}/ALBACORE_2.1.7_OUTPUT ${INPUT}/ALBACORE_2.3.4_OUTPUT ${INPUT}/GUPPY_2.3.1_OUTPUT ${INPUT}/GUPPY_3.0.3_OUTPUT; do
cd $d/epinano/minimap2/
head -1 per_site.var.csv.slided.onekmer.oneline.5mer.csv > per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv
awk -F"," "\$1 ~ /[^/$base/][^/$base/][/$base/][^/$base/][^/$base/]/" per_site.var.csv.slided.onekmer.oneline.5mer.csv >> per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv
mv per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv ${name}_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv #we change here the name according to the specific modification
cd $d/epinano/graphmap/
head -1 per_site.var.csv.slided.onekmer.oneline.5mer.csv > per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv
awk -F"," "\$1 ~ /[^/$base/][^/$base/][/$base/][^/$base/][^/$base/]/" per_site.var.csv.slided.onekmer.oneline.5mer.csv >> per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv
mv per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv ${name}_per_site.var.csv.slided.onekmer.oneline.5mer.filtered.csv
done
