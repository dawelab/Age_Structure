#!/bin/bash
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

