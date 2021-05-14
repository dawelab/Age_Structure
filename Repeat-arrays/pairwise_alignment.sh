#!/bin/bash
input_fa1="$1"
input_fa2="$2"

method=${3-water}
echo $method
#method="water"

module load SAMtools/1.9-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0

type=$(basename $input_fa1 |cut -f1 -d ".")
prefix=$(basename $input_fa1 |cut -f3 -d ".")
chr=$(basename $input_fa1 | cut -f2 -d ".")

line1=$(basename $input_fa1 | cut -f4 -d ".")
line2=$(basename $input_fa2 | cut -f4 -d ".")

wdir=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/pairwise_alignment/${type}/${chr}/${prefix}
alignment_output=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/pairwise_alignment/${type}/${chr}/${prefix}/${type}.${chr}.${prefix}.${line1}_${line2}.txt
plot_input=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/pairwise_alignment/${type}/${chr}/${prefix}/${type}.${chr}.${prefix}.${line1}_${line2}.labelled.stat.plot
plot_output=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/pairwise_alignment/${type}/${chr}/${prefix}/${type}.${chr}.${prefix}.${line1}_${line2}.png


runwater () {
    local run=$1
    water -asequence $x -bsequence $y -outfile align -gapopen 10 -gapextend 0.5 -outfile ${line1}-${xprefix}_${line2}-${yprefix}.water -brief Y
    aligned=$(cat ${line1}-${xprefix}_${line2}-${yprefix}.water | grep "Identity" | cut -f2 -d ":" | cut -f1 -d "/" | sed 's|[[:blank:]]||g')
    identity=$(cat ${line1}-${xprefix}_${line2}-${yprefix}.water | grep "Identity" | cut -f2 -d "(" | cut -f1 -d "%" | sed 's|[[:blank:]]||g')
    jaccard=$(echo "scale=2;100*$aligned/($xlen+$ylen-$aligned)" | bc )
    echo -e $xprefix_num"\t"$xlen"\t"$yprefix_num"\t"$ylen"\t"$aligned"\t"$jaccard"\t"$identity >> $alignment_output
    rm ${line1}-${xprefix}_${line2}-${yprefix}.water
}



mkdir -p $wdir
cd $wdir

if [[ $method == "water" ]]
then
  echo "EMBOSS WATER alignment"
  ml EMBOSS/6.6.0-GCC-8.3.0-Java-11
  ml bioawk/1.0-foss-2019b
for x in /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${prefix}/${line1}/*.fa
do
xprefix=$(basename $x |cut -f5,6 -d ".")
xprefix_num=$(basename $x |cut -f6 -d ".")
xlen=$(bioawk -c fastx '{ print length($seq) }' $x)
for y in /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${prefix}/${line2}/*.fa
do
yprefix=$(basename $y |cut -f5,6 -d ".")
yprefix_num=$(basename $y |cut -f6 -d ".")
ylen=$(bioawk -c fastx '{ print length($seq) }' $y)
runwater "$run" & done done
# water -asequence $x -bsequence $y -outfile align -gapopen 10 -gapextend 0.5 -outfile ${line1}-${xprefix}_${line2}-${yprefix}.water -brief Y
# aligned=$(cat ${line1}-${xprefix}_${line2}-${yprefix}.water | grep "Identity" | cut -f2 -d ":" | cut -f1 -d "/" | sed 's|[[:blank:]]||g')
# identity=$(cat ${line1}-${xprefix}_${line2}-${yprefix}.water | grep "Identity" | cut -f2 -d "(" | cut -f1 -d "%" | sed 's|[[:blank:]]||g')
# jaccard=$(echo "scale=2;100*$aligned/($xlen+$ylen-$aligned)" | bc )
# echo -e $xprefix_num"\t"$xlen"\t"$yprefix_num"\t"$ylen"\t"$aligned"\t"$jaccard"\t"$identity >> $alignment_output
# rm ${line1}-${xprefix}_${line2}-${yprefix}.water


# all by all sequence alignment result for plotting
cat $alignment_output | \
awk -v var1=$line1.${prefix} -v var3=$line2.${prefix} -v var2=$chr -v OFS='\t' '{print var2,var1,$1,var3,$3,$6}' | \
sort -k3,3n -k5,5n > $plot_input

else
  echo "BLAT"
  ml BLAT/3.5-GCC-8.3.0
  mkdir -p psl
  cd psl
  for x in /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${prefix}/${line1}/*.fa
  do
  xprefix=$(basename $x |cut -f5,6 -d ".")
  #xprefix_num=$(basename $x |cut -f6 -d ".")
  psl_output=${type}.${chr}.${prefix}.${line1}_${line2}.${xprefix}.psl
  blat -minIdentity=70 -maxGap=10 -minScore=0 -repMatch=2147483647 $x $input_fa2 $psl_output
  done
  cd ..
  cat psl/${type}.${chr}.${prefix}.${line1}_${line2}.*.psl > $alignment_output
  cat $alignment_output | grep $type | \
  awk -v OFS='\t' '{print $1,$14,$15,$10,$11}' | \
  #awk -v var1=$line1 -v var2=$line2 -v OFS='\t' '{print $1, var2, $2,$3, var1,$4,$5}' | \
  #awk -v OFS='\t' '{if($3>$6){print$1, $5,$6,$7,$2,$3,$4}else{print$0}}' | \
  awk -v OFS='\t' '{print$2, $3, $4, $5, $1, $1*100/($3+$5-$1)}' | \
  awk -v var1=$type. -v OFS='\t' 'BEGIN{a[1];a[3]}{for(x in a)gsub(var1,"",$x);{print}}' | \
  awk -v var1=$line1.${prefix} -v var3=$line2.${prefix} -v var2=$chr -v OFS='\t' '{print var2,var1,$1,var3,$3,$6}' | \
  sort -k3,3n -k5,5n -u > $plot_input
fi

# ml R
# Rscript ~/git/NAM_pancentromere/syn_repeat/plot_scatter.R \
# $wdir \
# $plot_input \
# $plot_output




#awk -v var1=$type -v OFS='\t' 'BEGIN{a[3];a[5]}{for(x in a)gsub(var1,"",$x);{print}}'
