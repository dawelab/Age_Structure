#!/bin/bash

chr="$1"
ml Anaconda3
# synaligned_file=/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/pangenome/NAM.$chr.pairwise.syn.bed
# cat /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/$chr/*.syn.aligned.bed > $synaligned_file
# python /home/jl03308/git/NAM_pancentromere/pangenome/pangenome_cal_shuf.py \
# -i $synaligned_file \
# -g /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/NAM.genome \
# -c $chr \
# -o /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/pangenome/$chr.pangenome.bed \
# -n 1000

synaligned_file=/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/pangenome/NAM.$chr.pairwise.syn.bed
cat /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/$chr/*.syn.aligned.bed > $synaligned_file
python /home/jl03308/git/NAM_pancentromere/pangenome/pangenome_cal_shuf_parallel.py \
-i $synaligned_file \
-g /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/NAM.genome \
-c $chr \
-o /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/pangenome/$chr.pangenome.bed \
-n 1000
