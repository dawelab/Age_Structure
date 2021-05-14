#!/bin/bash
type="$1"
chr="$2"
prefix="$3"

# type=knob180
# chr=chr9
# prefix=knob1

wdir=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/pairwise_alignment/$type/${chr}/$prefix
outdir=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating

distance_measure () {
    #local run=$1
    alignedlen=$(cat $file |awk -v var1=$query.$prefix '{if($6>=98.0&&$4==var1){print$0}}' | wc -l)
    if [ "$alignedlen" -gt 0 ]; then
  cat $file | \
  awk -v var1=$query.$prefix '{if($6>=98.0&&$4==var1){print$0}}' | \
  sort -k3,3n -k5,5n | awk -F"\t" '!seen[$3]++' | \
  sort -k5,5n -k3,3n | awk -F"\t" '!seen[$5]++'| \
  awk -v var1=$ref_line -v var2=$query -v var3=$chr -v var4=$prefix -v var5=$ref_seq_len -v var6=$query_seq_len -v var7=$type -v OFS='\t' '{sum+=$6} END {print var3, var4, var7, var1, var2, var5, var6, sum, NR, sum/NR}' >> ${outdir}/NAM.$type.$chr.$prefix.pairwise_distance.txt
  else
    echo -e $chr"\t"$prefix"\t"$type"\t"$ref_line"\t"$query"\t"$ref_seq_len"\t"$query_seq_len"\t0\t0\t0" >> ${outdir}/NAM.$type.$chr.$prefix.pairwise_distance.txt
  fi
}

for file in ${wdir}/${type}.${chr}.*_*.labelled.stat.plot
do
#file=knob180.chr9.knob1.Tx303.ref.labelled.stat.plot
ref_line=$(basename $file |cut -f4 -d "." |cut -f1 -d "_")
query=$(basename $file |cut -f4 -d "." |cut -f2 -d "_")
refseq_bed_file=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${ref_line}.bed
ref_seq_len=$(cat $refseq_bed_file | tail -1 |cut -f4 |cut -f2 -d ".")
queryseq_bed_file=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${query}.bed
query_seq_len=$(cat $queryseq_bed_file | tail -1 |cut -f4 |cut -f2 -d ".")
distance_measure #"$run" &
done

ml R
Rscript ~/git/NAM_pancentromere/syn_repeat/knob_cluster.R  \
$outdir \
${outdir}/NAM.${type}.${chr}.${prefix}.pairwise_distance.txt \
NAM.${type}.${chr}.${prefix}
