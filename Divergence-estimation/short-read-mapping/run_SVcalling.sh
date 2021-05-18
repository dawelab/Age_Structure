#!/bin/bash

bam="$1"
ref_genome="$2"

prefix=$(basename $bam | cut -f1,2 -d ".")
wdir=$(dirname $bam)

bin=1000

cd $wdir
sh /home/jl03308/git/NAM_pancentromere/teosinte/mpileup.sh $bam ${ref_genome}
sh /home/jl03308/git/NAM_pancentromere/teosinte/deeptools.sh $bam $bin
