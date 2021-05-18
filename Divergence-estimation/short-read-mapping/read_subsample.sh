#!/bin/bash
PE1="$1"
PE2="$2"
num="$3"

module load seqtk/1.3-GCC-8.3.0
wdir=$(dirname $PE1)
prefix=$(basename $PE1 | cut -f1 -d ".")
prefix2=$(basename $PE2 | cut -f1 -d ".")

seqtk sample -2 -s 11 $PE1 $num |gzip > ${wdir}/${prefix}.${num}.fastq.gz
seqtk sample -2 -s 11 $PE2 $num |gzip > ${wdir}/${prefix2}.${num}.fastq.gz
