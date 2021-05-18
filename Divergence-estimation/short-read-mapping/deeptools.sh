#!/bin/bash

bam="$1"
bin="$2"

#bam=/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/TIL07/TIL07.B73.q20.sorted.bam
#bin=1000

wdir=$(dirname $bam)
prefix=$(basename $bam |cut -f1,2 -d ".")
cd $wdir

ml deepTools/3.3.1-intel-2019b-Python-3.7.4
#bamCoverage -b $bam -p 12 -bs $bin -o ${prefix}.bedgraph -of bedgraph --minMappingQuality 20 --ignoreDuplicates
bamCoverage -b $bam -p 12 -bs $bin -o ${prefix}.RPKM.bedgraph -of bedgraph --minMappingQuality 20 --ignoreDuplicates --normalizeUsing RPKM
