#!/bin/bash
RAW=$1                # directory with folders with single fast5s, or directory with multifast5s
OUTPUT=$2             # name of the directory that will be created with all the outputs
multifast5=$3         # boolean, true or false

source /users/enovoa/joramirez/.bashrc #We source in José's bashrc file to have some of my software available


if [[ ! -e ${OUTPUT} ]]
then
 mkdir $OUTPUT
 mkdir ${OUTPUT}/GUPPY_3.0.3_OUTPUT
 mkdir ${OUTPUT}/GUPPY_2.3.1_OUTPUT
 mkdir ${OUTPUT}/ALBACORE_2.3.4_OUTPUT
 mkdir ${OUTPUT}/ALBACORE_2.1.7_OUTPUT
else
 echo "[ WARNING ] Found existing directories"
fi



if $multifast5 ; then
 mkdir -p ${OUTPUT}/singleRAW
 /users/enovoa/joramirez/scripts/launch.ont-basecalling_guppy_multifast5.sh $RAW ${OUTPUT}/GUPPY_3.0.3_OUTPUT
 /users/enovoa/joramirez/scripts/launch.ont-basecalling_guppy_multifast5_2.3.1.sh $RAW ${OUTPUT}/GUPPY_2.3.1_OUTPUT
 #multi_to_single_fast5 -i $RAW -s ${OUTPUT}/singleRAW/ -t 8
 pip3 install /users/enovoa/joramirez/software/ont_albacore-2.1.7-cp36-cp36m-manylinux1_x86_64.whl --user #installing ALBACORE 2.1.7 in your user
 /users/enovoa/joramirez/scripts/launch.ont-basecalling_untared.sh ${OUTPUT}/singleRAW/ ${OUTPUT}/ALBACORE_2.1.7_OUTPUT
 pip3 install /users/enovoa/joramirez/software/ont_albacore-2.3.4-cp36-cp36m-manylinux1_x86_64.whl --user #installing ALBACORE 2.3.4 in your user
 /users/enovoa/joramirez/scripts/launch.ont-basecalling_untared.sh ${OUTPUT}/singleRAW/ ${OUTPUT}/ALBACORE_2.3.4_OUTPUT
else
 /users/enovoa/joramirez/scripts/launch.ont-basecalling_guppy.sh $RAW ${OUTPUT}/GUPPY_3.0.3_OUTPUT
 /users/enovoa/joramirez/scripts/launch.ont-basecalling_guppy_2.3.1.sh $RAW ${OUTPUT}/GUPPY_2.3.1_OUTPUT
 pip3 install /users/enovoa/joramirez/software/ont_albacore-2.1.7-cp36-cp36m-manylinux1_x86_64.whl --user
 /users/enovoa/joramirez/scripts/launch.ont-basecalling.sh $RAW ${OUTPUT}/ALBACORE_2.1.7_OUTPUT
 pip3 install /users/enovoa/joramirez/software/ont_albacore-2.3.4-cp36-cp36m-manylinux1_x86_64.whl --user
 /users/enovoa/joramirez/scripts/launch.ont-basecalling.sh $RAW ${OUTPUT}/ALBACORE_2.3.4_OUTPUT
fi





