## Tools used for SV calling and alignment plotting
- Alignment
    * bedtools (2020.02)
    * BLAT (v3.5)

- Plot
    * qgraph

## Identify syntenic arrays
Format example `chr3.knob1.bed`:
```
B97	186255901	189405832
CML247	191211695	192938110
CML277	192050552	193091267
CML333	192358172	193129332
HP301	189584799	192175719
M37W	192619570	193229744
MS71	189982971	192390352
Mo18W	191406204	193490188
NC350	190160515	194617525
```

## Extract repeat sequences
Run `repeat_extract.sh` to extract, and index monomers in each array
```
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
```

## Pairwise alignment
Run `pairwise_alignment.sh` for each pairwise alignment of two arrays
```
syn_repeat_file="$1"

chr=$(basename $syn_repeat_file | cut -f1 -d ".")
prefix=$(basename $syn_repeat_file | cut -f2 -d ".")
repeat_type=$(basename $syn_repeat_file | cut -f2 -d "." |sed 's/[0-9]*//g')

if [[ $repeat_type == "knob" ]]
then
  types=(knob180 TR-1)
else
  types=($repeat_type)
fi

for type in "${types[@]}"
do
  echo $type
#pairwise-alignment of knobs
all_line=($(cat $syn_repeat_file |cut -f1 |sort |uniq))
len=${#all_line[@]}
## Use bash for loop
for ((i = 0; i < len; i++)); do
  for ((j = i+1; j < len; j++)); do
    line1=${all_line[$i]}
    line2=${all_line[$j]}
    #echo "$line1,$line2"
    input_fa1=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${line1}.fa
    input_fa2=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${line2}.fa
    sh pairwise_alignment.sh $input_fa1 $input_fa2
  done
done
done
```

## Plot

```
type="$1"
chr="$2"
prefix="$3"

# type=knob180
# chr=chr5
# prefix=knob1

ml R
wdir=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/pairwise_alignment/${type}/${chr}/${prefix}
cd $wdir
repeat_array_file=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/arrays/${chr}.${prefix}.bed
NAM_lines=$(cat $repeat_array_file | cut -f1 | sort |uniq)
for line in $NAM_lines
do
  plot_input=${type}.${chr}.${prefix}.${line}.ref.labelled.stat.plot
  if test -f "$plot_input"; then
    echo "$plot_input exists."
  else
  cat ${type}.${chr}.${prefix}.${line}_*.labelled.stat.plot >> $plot_input
  cat ${type}.${chr}.${prefix}.*_${line}.labelled.stat.plot | awk -v OFS='\t' '{print $1,$4,$5,$2,$3,$6}' >> $plot_input
fi
  plot_output=${type}.${chr}.${prefix}.${line}.ref.png
  Rscript ~/git/NAM_pancentromere/syn_repeat/plot_scatter.R \
  $wdir \
  $plot_input \
  $plot_output
done
```

