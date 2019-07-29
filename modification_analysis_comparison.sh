base=$1
unm=$2
mod=$3
OUTPUT=$4
model=$5

echo $model

#filtering needed for building epinano model
for i in $unm $mod; do 
head -1 $i > ${i%.*}.filtered.csv
awk -F"," "\$1 ~ /[^/$base/][^/$base/][/$base/][^/$base/][^/$base/]/" $i >> ${i%.*}.filtered.csv
done

if [[ $model == true ]]; then

echo Building epinano model...

cut -d "," -f 5-14,20-24 ${unm} > test
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

cut -d "," -f 5-14,20-24 ${mod} > test
tail -n +2 test > test2
sed 's/$/,MOD/' test2 > test3
head -$t test3 >> training 
tail -$p test3 >> predicting
rm test test2 test3 

SVM.py -a -t training -p predicting -cl 1-15 -mc 16 -o $OUTPUT/svm_$base

elif [[ $model == false ]]; then

echo Skipping epinano model

else
echo No model option provided; fi




echo Plotting...

if [ -z "$6" ]; then
	echo Using "UNM,MOD" as dataset names
	Rscript scripts/kmer.R -u ${unm%.*}.filtered.csv -m ${mod%.*}.filtered.csv -o $OUTPUT
else
	Rscript scripts/kmer.R -u ${unm%.*}.filtered.csv -m ${mod%.*}.filtered.csv -o $OUTPUT -n $6
fi


