#!/bin/bash

syn_repeat_file="$1"

chr=$(basename $syn_repeat_file | cut -f1 -d ".")
prefix=$(basename $syn_repeat_file | cut -f2 -d ".")
repeat_type=$(basename $syn_repeat_file | cut -f2 -d "." |sed 's/[0-9]*//g')

#extract knob sequences
if [[ $repeat_type == "knob" ]]
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
