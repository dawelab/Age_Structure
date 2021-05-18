#!/bin/bash

bam="$1"
ref_genome="$2"
#ref_genome=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/B73.PLATINUM.pseudomolecules-v1.fasta
#bam=/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/TIL07/TIL07.B73.q20.sorted.bam

wdir=$(dirname $bam)
prefix=$(basename $bam |cut -f1,2 -d ".")
cd $wdir

ml SAMtools/1.9-GCC-8.3.0
ml BCFtools/1.6-foss-2019b

bcftools mpileup -Ou -C50 -f $ref_genome $bam | \
bcftools call -Ou -mv --threads 12 > ${prefix}.raw.bcf
dp=$(samtools depth $bam |  awk '{sum+=$3} END { print sum/NR}')
dp_cutoff_up=$(echo "scale=1;$dp*4" | bc)
dp_cutoff_low=$(echo "scale=1;$dp/4" | bc)

bcftools filter -Oz -e "%QUAL<20 || DP>$dp_cutoff_up || DP<$dp_cutoff_low" ${prefix}.raw.bcf > ${prefix}.flt.vcf.gz

bcftools view -v snps -g hom ${prefix}.flt.vcf.gz -Ou | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\t%QUAL\n' | \
awk -v OFS='\t' '{print $1,$2,$2+1,$3,$4,$5}'> ${prefix}.flt.snps.bed
