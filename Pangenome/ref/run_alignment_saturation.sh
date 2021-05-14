#!/bin/bash
chr="$1"

ml Anaconda3
wdir=/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/$chr
cat $wdir/B73_*.syn.aligned.bed > $wdir/NAM_toB73.syn.aligned.bed
python /home/jl03308/git/NAM_pancentromere/pangenome_ref/alignment_saturation.py \
-i $wdir/NAM_toB73.syn.aligned.bed \
-o /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/$chr.alignment_saturation.bed \
-n 1000
rm $wdir/NAM_toB73.syn.aligned.bed
