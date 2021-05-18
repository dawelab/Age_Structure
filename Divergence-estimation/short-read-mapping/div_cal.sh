#!/bin/bash

effectivecov_file="$1"
snp_file="$2"
ref="$3"

prefix=$(basename $effectivecov_file | cut -f1,2 -d ".")
outdir=$(dirname $effectivecov_file)

ml BEDTools/2.29.2-GCC-8.3.0
#ref=B73
window_size=20000
window_file=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.${window_size}.windows
wdir=/scratch/jl03308/NAM_pancentromere/teosinte_divergence

#snp_file=/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/PC_O59_ID2/PC_O59_ID2.q20.flt.snps.bed
#effectivecov_file=/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/PC_O59_ID2/PC_O59_ID2.q20.effective.genomecov


# cd $wdir
# # gene and UMR file
# for chr in chr{1..10}
# do
#   gene_UMR_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency/${ref}.${chr}.UMR_gene.bed
#   cat $gene_UMR_file | awk -v var1=$chr '{print var1"\t"$2"\t"$3}' >> $wdir/$ref.UMR_gene.bed
# done

# identify snps not in genes or UMR
bedtools subtract -a $snp_file -b $wdir/${ref}.UMR_gene.bed | \
bedtools sort -i - > ${outdir}/${prefix}.flt.snps.nonUMR_gene.bed

# remove effective coverage in gene and umr regions
cat $effectivecov_file| \
bedtools merge -i - | \
bedtools subtract -a - -b $wdir/${ref}.UMR_gene.bed | \
bedtools sort -i - | \
# calculate effective coverage in window (no gene/UMR)
bedtools intersect -a $window_file -b - -wao| \
bedtools groupby -i - -g 1,2,3 -c 7 -o sum | \
# calculate number of snps in window (no gene/UMR)
bedtools intersect -a - -b ${outdir}/${prefix}.flt.snps.nonUMR_gene.bed -wao | \
bedtools groupby -i - -g 1,2,3,4 -c 11 -o sum | \
awk -v OFS='\t' '{if($4>0){print $0, $5/$4/2/3.3*100000000}else{print $0, "NA"}}' > ${outdir}/${prefix}.20Kb.div.txt
