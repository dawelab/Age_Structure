ml BEDTools/2.29.2-GCC-8.3.0
ref=B73
NAMline=(B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)

#identify all syntenic snps
for chr in chr{1..10}
do
  echo $chr
  out_wdir=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/$chr
  mkdir -p $out_wdir
for line in ${NAMline[@]}
do
  echo $line
  cat /scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window/${line}_${chr}.mapped-to.${ref}_${chr}-*_*-*_*.syn.snp.vcf | \
  grep -v "#" | \
  sed 's|=|\t|g' | sed 's|;|\t|g' |
  awk -v var1=$line -v var2=$chr -v OFS='\t' -v var3=$ref '{print var3,$2,$2+1,var1,$11,$11+1,var2,$4,$5}' | \
  bedtools sort -i - |uniq > ${out_wdir}/${line}.${chr}.syn.snp.bed
done
done

#map all syntenic snps to intergenic locations on each aligned blocks

for chr in chr{1..10}
do
  for line in ${NAMline[@]}
  do
    echo $line
    syn_snp_file=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/$chr/${line}.${chr}.syn.snp.bed
    gene_UMR_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency/${ref}.${chr}.UMR_gene.bed
    syn_alignment_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/$chr/${ref}_${line}.syn.aligned.bed
  # remove syn snp from gene and UMR
  bedtools intersect -a $syn_snp_file -b $gene_UMR_file -v | \
  # map syn snps in intergenic space to syn aligned blocks
  bedtools intersect -a $syn_alignment_file -b - -c| \
  # calculate intergenic space in each aligned blocks
  bedtools intersect -a - -b $gene_UMR_file -wao | \
  bedtools sort -i - | \
  # print intergenic syn snps, and intergenic length in each aligned blocks
  # ref, start, end, query, start,end, intergenic snp number, gene+UMR length
  bedtools groupby -i - -g 1,2,3,4,5,6,14 -c 18 -o sum | \
  # ref, start, end, query, start,end, intergenic snp number, intergenic length
  awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$3-$2-$8}' | \
  awk -v OFS='\t' -v var1=$chr '{if($8>0){print $0,$7/$8/2/3.3 *100000000,var1}}' >> /scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_div.bed
done
done

#map all syntenic snps to intergenic locations on each window
window_size=20000
window_file=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.${window_size}.windows
for chr in chr{1..10}
do
  echo $chr
for line in ${NAMline[@]}
do
  echo $line
  syn_snp_file=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/$chr/${line}.${chr}.syn.snp.bed
  gene_UMR_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency/${ref}.${chr}.UMR_gene.bed
  syn_alignment_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/$chr/${ref}_${line}.syn.aligned.bed
  #calculate non-gene/UMR regions in each window
cat $window_file | \
awk -v var1=$chr -v var2=$ref -v OFS='\t' '{if($1==var1){print var2,$2,$3}}' | \
bedtools intersect -a - -b $syn_alignment_file -wao | \
bedtools sort -i - | \
bedtools groupby -i - -g 1,2,3,7 -c 17 -o sum | \
bedtools intersect -a - -b $gene_UMR_file -wao | \
bedtools sort -i - | \
# print intergenic syn snps, and intergenic length in each aligned blocks
# ref, start, end, query, start,end, intergenic snp number, gene+UMR length
bedtools groupby -i - -g 1,2,3,4,5 -c 9 -o sum | \
awk -v var1=$line -v var2=$chr -v OFS='\t' '{if($5>0){print $1,$2,$3,var1,$5-$6,var2}else{print $1,$2,$3,var1,"NA",var2}}' >> /scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_${window_size}_aligned.bed
done
done

#map snps on non-gene/UMR regions to each window
for chr in chr{1..10}
do
  echo $chr
for line in ${NAMline[@]}
do
  syn_snp_file=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/$chr/${line}.${chr}.syn.snp.bed
  gene_UMR_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency/${ref}.${chr}.UMR_gene.bed
  syn_alignment_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/$chr/${ref}_${line}.syn.aligned.bed
  cat /scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_${window_size}_aligned.bed | \
  awk -v var1=$line -v var2=$chr -v OFS='\t' '{if($4==var1&&$6==var2){print $1,$2,$3,$4,$5,$6}}' > /scratch/jl03308/NAM_pancentromere/divergence/intergenic/tmp/$chr.$line.intergenic_${window_size}_aligned.bed
  bedtools intersect -a $syn_snp_file -b $gene_UMR_file -v | \
  # map syn snps in intergenic space to syn aligned blocks
  bedtools intersect -a /scratch/jl03308/NAM_pancentromere/divergence/intergenic/tmp/$chr.$line.intergenic_${window_size}_aligned.bed -b - -c| \
  awk -v OFS='\t' '{if($5!="NA"&&$5>0){print $0,$7/$5/2/3.3 *100000000}}' >> /scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_${window_size}_div.bed
done
done


#map all syntenic snps to intergenic locations according to allele frequency

allele_frequency_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/alignment_frequency_allchrs.bed
for num in {1..25}
do
for chr in chr{1..10}
do
  echo $chr
  gene_UMR_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency/${ref}.${chr}.UMR_gene.bed
cat $allele_frequency_file | \
awk -v var1=$chr -v var2=$ref -v var3=$num -v OFS='\t' '{if($28==var1&&$29==var3){print var2,$1,$2}}' | \
bedtools sort -i - | \
bedtools merge -i -| \
bedtools intersect -a - -b $gene_UMR_file -wao | \
bedtools sort -i - | \
# print intergenic syn snps, and intergenic length in each aligned blocks
# ref, start, end, query, start,end, intergenic snp number, gene+UMR length
bedtools groupby -i - -g 1,2,3 -c 7 -o sum | \
awk -v var2=$chr -v OFS='\t' '{print var2,$1,$2,$3,$4,$3-$2-$4}' >> /scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_${num}.bed
done
done

for num in {1..25}
do
for chr in chr{1..10}
do
  echo $chr
for line in ${NAMline[@]}
do
  syn_snp_file=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/$chr/${line}.${chr}.syn.snp.bed
  gene_UMR_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency/${ref}.${chr}.UMR_gene.bed
  cat /scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_${num}.bed| \
  awk -v var2=$chr -v OFS='\t' '{if($1==var2){print $2,$3,$4,$5,$6}}' > /scratch/jl03308/NAM_pancentromere/divergence/intergenic/tmp/NAM.intergenic_${num}.$chr.bed
  bedtools intersect -a $syn_snp_file -b $gene_UMR_file -v | \
  # map syn snps in intergenic space to syn aligned blocks
  bedtools intersect -a /scratch/jl03308/NAM_pancentromere/divergence/intergenic/tmp/NAM.intergenic_${num}.$chr.bed -b - -c| \
  awk -v OFS='\t' -v var1=$chr -v var2=$num '{if($5!="NA"&&$5>50){print $0,$6/$5/2/3.3 *100000000,var1,var2}}' >> /scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_pangenome_div.bed
done
done
done
