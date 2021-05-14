#!/bin/bash

syn_repeat_file="$1"

chr=$(basename $syn_repeat_file | cut -f1 -d ".")
prefix=$(basename $syn_repeat_file | cut -f2 -d ".")
repeat_type=$(basename $syn_repeat_file | cut -f2 -d "." |sed 's/[0-9]*//g')

#extract knob sequences
if [ $repeat_type == "knob" ]
then
  types=(knob180 TR-1)
else
  types=($repeat_type)
fi

for type in "${types[@]}"
do
  echo $type
while IFS= read line
do
  if [ ! -z "$line" ]; then
nam_line=$(echo $line |cut -f1 -d " ")
start=$(echo $line |cut -f2 -d " ")
end=$(echo $line |cut -f3 -d " ")
echo -e $nam_line"\t"$chr"\t"$type"\t"$prefix"\t"$start"\t"$end
sh /home/jl03308/git/NAM_pancentromere/syn_repeat/repeat_extract.sh \
$nam_line \
$chr \
$start \
$end \
$type \
$prefix
fi
done < $syn_repeat_file
done

#self-alignment of knobs
while IFS= read line
do
  if [ ! -z "$line" ]; then
nam_line=$(echo $line |cut -f1 -d " ")
input_fa=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${nam_line}.fa
echo -e $input_fa
sh /home/jl03308/git/NAM_pancentromere/syn_repeat/self_alignment.sh $input_fa
fi
done < $syn_repeat_file

#pairwise-alignment of knobs
all_line=($(cat $syn_repeat_file |cut -f1 |sort |uniq))
len=${#all_line[@]}
## Use bash for loop
for ((i = 0; i < len; i++)); do
  for ((j = i+1; j < len; j++)); do
    line1=${all_line[$i]}
    line2=${all_line[$j]}
    echo "$line1,$line2"
    input_fa1=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${line1}.fa
    input_fa2=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${line2}.fa
    sh /home/jl03308/git/NAM_pancentromere/syn_repeat/pairwise_alignment.sh $input_fa1 $input_fa2
  done
done
