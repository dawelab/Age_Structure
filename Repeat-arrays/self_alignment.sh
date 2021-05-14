#!/bin/bash
input_fa="$1"

module load SAMtools/1.9-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0
ml EMBOSS/6.6.0-GCC-8.3.0-Java-11
ml bioawk/1.0-foss-2019b

type=$(basename $input_fa |cut -f1 -d ".")
line=$(basename $input_fa | cut -f4 -d ".")
chr=$(basename $input_fa | cut -f2 -d ".")
prefix=$(basename $input_fa |cut -f3 -d ".")


wdir=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/self_alignment/${type}/${chr}/${line}
alignment_output=${type}.${chr}.${prefix}.${line}.txt
plot_input=${type}.${chr}.${prefix}.${line}.labelled.stat.plot
plot_output=${type}.${chr}.${prefix}.${line}.png

mkdir -p $wdir
cd $wdir

for x in /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${prefix}/${line}/*.fa
do
xprefix=$(basename $x |cut -f5,6 -d ".")
xprefix_num=$(basename $x |cut -f6 -d ".")
xlen=$(bioawk -c fastx '{ print length($seq) }' $x)
for y in /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${prefix}/${line}/*.fa
do
yprefix=$(basename $y |cut -f5,6 -d ".")
yprefix_num=$(basename $y |cut -f6 -d ".")
ylen=$(bioawk -c fastx '{ print length($seq) }' $y)
if [ "$xprefix_num" -gt "$yprefix_num" ]
then
water -asequence $x -bsequence $y -outfile align -gapopen 10 -gapextend 0.5 -outfile ${xprefix}_${yprefix}.water -brief Y
aligned=$(cat ${xprefix}_${yprefix}.water | grep "Identity" | cut -f2 -d ":" | cut -f1 -d "/" | sed 's|[[:blank:]]||g')
identity=$(cat ${xprefix}_${yprefix}.water | grep "Identity" | cut -f2 -d "(" | cut -f1 -d "%" | sed 's|[[:blank:]]||g')
jaccard=$(echo "scale=2;100*$aligned/($xlen+$ylen-$aligned)" | bc )
echo -e $xprefix_num"\t"$xlen"\t"$yprefix_num"\t"$ylen"\t"$aligned"\t"$jaccard"\t"$identity >> $alignment_output
rm ${xprefix}_${yprefix}.water
fi
done
done

# all by all sequence alignment result for plotting
cat $alignment_output | \
awk -v var1=$line.${prefix} -v var2=$chr -v OFS='\t' '{print var2,var1,$1,var1,$3,$6}' | \
sort -k3,3n -k5,5n > $plot_input

ml R
Rscript ~/git/NAM_pancentromere/syn_repeat/plot_scatter.R \
$wdir \
$plot_input \
$plot_output


#awk -v var1=$type -v OFS='\t' 'BEGIN{a[3];a[5]}{for(x in a)gsub(var1,"",$x);{print}}'
