#!/bin/bash
chr="$1"
ref="$2"
ml BEDTools/2.29.2-GCC-8.3.0

NAMline=(B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
ref_lc=$(echo "$ref" | awk '{print tolower($0)}')
wdir=/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency
#chr=chr8
#ref=B73
chrlen=$(cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/$ref.*pseudomolecules-v*.fasta.fai | awk -v var1=$chr '{if($1==var1){print $2}}')
cd $wdir
cat ${chr}.allele_frequency.bed | \
awk -v var1=$ref -v OFS='\t' '{print var1,$0,$2-$1+1}' | \
sort -k1,1 -k4,4n | \
bedtools groupby -i - -g 1,4 -o sum -c 5 > ${chr}.allele_frequency.totallen.bed

#genes
cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/${ref_lc}/zea_mays${ref_lc}_core_3_87_1.gff | \
awk -v var1=$chr -v var2=$ref -v OFS='\t' '{if($3=="gene"&&$1==var1){print var2,$4,$5}}' | \
bedtools sort -i - | \
bedtools merge -i - > ${ref}.${chr}.gene.bed
#UMRs
cat /scratch/jl03308/NAM_pancentromere/methylation/UMRs/meth_${ref}.ref_${ref}.UMR.bed | \
awk '{print "chr"$0}' | \
awk -v var1=$chr -v var2=$ref -v OFS='\t' '{if($1==var1){print var2,$2,$3}}' | \
bedtools sort -i - | \
bedtools merge -i - > ${ref}.${chr}.UMR.bed

cat ${ref}.${chr}.gene.bed ${ref}.${chr}.UMR.bed |bedtools sort -i - | bedtools merge -i - > ${ref}.${chr}.UMR_gene.bed

totalgenelen=$(cat ${ref}.${chr}.gene.bed | awk '{sum+=$3-$2+1;} END{print sum;}')
totalUMRlen=$(cat ${ref}.${chr}.UMR.bed | awk '{sum+=$3-$2+1;} END{print sum;}')
totalgeneUMRlen=$(cat ${ref}.${chr}.UMR_gene.bed | awk '{sum+=$3-$2+1;} END{print sum;}')

cat ${chr}.allele_frequency.bed | \
awk -v var1=$ref -v OFS='\t' '{print var1,$0}' | \
bedtools intersect -a - -b ${ref}.${chr}.gene.bed -wao | \
sort -k1,1 -k4,4n | \
bedtools groupby -i - -g 1,4 -o sum -c 8 > ${chr}.allele_frequency.genelen.bed

cat ${chr}.allele_frequency.bed | \
awk -v var1=$ref -v OFS='\t' '{print var1,$0}' | \
bedtools intersect -a - -b ${ref}.${chr}.UMR.bed -wao | \
sort -k1,1 -k4,4n | \
bedtools groupby -i - -g 1,4 -o sum -c 8 > ${chr}.allele_frequency.UMRlen.bed

cat ${chr}.allele_frequency.bed | \
awk -v var1=$ref -v OFS='\t' '{print var1,$0}' | \
bedtools intersect -a - -b ${ref}.${chr}.UMR_gene.bed -wao | \
sort -k1,1 -k4,4n | \
bedtools groupby -i - -g 1,4 -o sum -c 8 > ${chr}.allele_frequency.UMR_genelen.bed

stat_file=/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency/$ref.${chr}.stat
# if test -f "$stat_file"; then
#   rm $stat_file
# else
paste ${chr}.allele_frequency.totallen.bed ${chr}.allele_frequency.genelen.bed ${chr}.allele_frequency.UMRlen.bed ${chr}.allele_frequency.UMR_genelen.bed| \
awk -v OFS='\t' -v var1=$totalgenelen -v var2=$totalUMRlen -v var3=$chrlen -v var4=$totalgeneUMRlen -v var5=$chr \
'{print $1,$2,$3,$3/var3*100,$6,$6/$3*100,$6/var1*100,$9,$9/$3*100,$9/var2*100,$12,$12/$3*100,$12/var4*100, var5}' >> $stat_file
#fi
