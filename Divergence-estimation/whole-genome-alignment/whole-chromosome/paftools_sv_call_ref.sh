#!/bin/bash
chr="$1"
ref="$2"
synpaf_file="$3"

ref_genome=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.fasta
prefix=$(basename $synpaf_file |cut -f1-4 -d ".")
#chr=chr8
#genome_query=B97
#module load BEDTools/2.29.2-GCC-8.3.0
ml minimap2
cat $synpaf_file | sed 's|B73_||g' | sort -k8,8n | /home/jl03308/bin/k8 /home/jl03308/minimap2/misc/paftools.js call -L50 -q0 -l50 -f $ref_genome - > ${prefix}.vcf
