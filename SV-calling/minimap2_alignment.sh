#!/bin/bash
ref="$1"
query="$2"
out_dir="$3"

ref_genome=$(basename $ref| cut -f1 -d ".")
query_genome=$(basename $query| cut -f1 -d ".")
chr=$(basename $ref| cut -f2 -d ".")

module load minimap2
minimap2 -c -cx asm5 --no-kalloc --print-qname --cs=long -t 12 $ref $query > ${out_dir}/${query_genome}_${chr}.mapped-to.${ref_genome}_${chr}.paf
