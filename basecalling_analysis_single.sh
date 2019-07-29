#This script outputs the number of base-called reads from a given output fastq.gz files 
FASTQ_FILE=$1

n=$( awk '{if(NR%4==1) print $1}' ${FASTQ_FILE} | sed -e "s/^@//" | wc -l )
echo There are $n reads in the oinput file


#Add plotting read length and qualities too
