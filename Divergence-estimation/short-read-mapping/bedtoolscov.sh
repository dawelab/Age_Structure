#!/bin/bash

bam="$1"
#bam=/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/TIL07/TIL07.B73.q20.sorted.bam
#bin=1000

wdir=$(dirname $bam)
prefix=$(basename $bam |cut -f1,2 -d ".")
cd $wdir

module load BEDTools/2.29.2-GCC-8.3.0
#bamCoverage -b $bam -p 12 -bs $bin -o ${prefix}.bedgraph -of bedgraph --minMappingQuality 20 --ignoreDuplicates
#bamCoverage -b $bam -p 12 -bs $bin -o ${prefix}.RPKM.bedgraph -of bedgraph --minMappingQuality 20 --ignoreDuplicates --normalizeUsing RPKM
bedtools genomecov -ibam $bam -bga > $prefix.genomecov


ml SAMtools/1.9-GCC-8.3.0

dp=$(samtools depth $bam |  awk '{sum+=$3} END { print sum/NR}')
dp_cutoff_up=$(echo "scale=1;$dp*4" | bc)
dp_cutoff_low=$(echo "scale=1;$dp/4" | bc)

cat $prefix.genomecov | \
awk -v var1=$dp_cutoff_up -v var2=$dp_cutoff_low '{if($4<var1&&$4>var2){print$0}}' > $prefix.effective.genomecov
