ml BEDTools/2.29.2-GCC-8.3.0

nam=(B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
outdir=/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/SV_stats

for line in "${nam[@]}"
do
#cat /scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window/${line}_${chr}.mapped-to.${ref}_${chr}-*_*-*_*.syn.snp.vcf | sort -k2,2n > ${vcf_dir}/${line}/${line}.${chr}.syn.vcf
cd /lustre2/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/$line
for chr in chr{1..10}
do
  illumina_snp_num=$(zcat $line.${chr}.q20.flt.vcf.gz | grep -v "#"|wc -l)
  syn_snp_num=$(zcat $line.${chr}.syn.header.vcf.gz | grep -v "#"|wc -l)
  overlap=$(bedtools intersect -a $line.${chr}.q20.flt.vcf.gz -b $line.${chr}.syn.header.vcf.gz | wc -l)
  echo -e $line"\t"$chr"\t"$illumina_snp_num"\t"$syn_snp_num"\t"$overlap >> $outdir/NAM.syn-illumina_snp.txt
done
done


for line in "${nam[@]}"
do
#cat /scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window/${line}_${chr}.mapped-to.${ref}_${chr}-*_*-*_*.syn.snp.vcf | sort -k2,2n > ${vcf_dir}/${line}/${line}.${chr}.syn.vcf
cd /lustre2/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/$line
for chr in chr{1..10}
do
  illumina_mapping=$(zcat $line.mask.bed.gz | awk -v var1=$chr '{if($1==var1){print$0}}' | sort -k2,2n |bedtools merge -i - | awk '{sum+=$3-$2;} END{print sum;}')
  syn_mapping=$(cat /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/$chr/${line}.mask.bed | awk -v var1=$chr '{if($1==var1){print$0}}' | sort -k2,2n |bedtools merge -i - | awk '{sum+=$3-$2;} END{print sum;}')
  overlap=$(bedtools intersect -a /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/$chr/${line}.mask.bed -b $line.mask.bed.gz | awk -v var1=$chr '{if($1==var1){print$0}}' | sort -k2,2n |bedtools merge -i - | awk '{sum+=$3-$2;} END{print sum;}')
  echo -e $line"\t"$chr"\t"$illumina_mapping"\t"$syn_mapping"\t"$overlap >> $outdir/NAM.syn-illumina_mapping.txt
done
done


for line in "${nam[@]}"
do
cat $outdir/NAM.syn-illumina_snp.txt | awk -v var1=$line -v OFS='\t' '{if($1==var1){print$0}}' | \
awk -v var1=$line -v OFS='\t' '{sum1+=$3;sum2+=$4;sum3+=$5} END{print var1,sum1,sum2,sum3}' >> $outdir/NAM.syn-illumina_snp.wholegenome.txt
cat $outdir/NAM.syn-illumina_mapping.txt | awk -v var1=$line -v OFS='\t' '{if($1==var1){print$0}}' | \
awk -v var1=$line -v OFS='\t' '{sum1+=$3;sum2+=$4;sum3+=$5} END{print var1,sum1,sum2,sum3}' >> $outdir/NAM.syn-illumina_mapping.wholegenome.txt
done
