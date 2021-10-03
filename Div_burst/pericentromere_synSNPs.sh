#!/bin/bash
ref="$1"
chr="$2"
start_coord="$3"
end_coord="$4"

module load BEDTools/2.30.0-GCC-8.3.0

cd /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}
'''generate pericentromere bed file'''
echo -e $ref"\t"$start_coord"\t"$end_coord > ${ref}.pericentromere.bed
div_files=$(ls ${ref}_*.aligned.div.txt |grep -v "pericentromere")
#div_files=$(ls ${ref}_*.syn.aligned.bed |grep -v "pericentromere")

'''identify syntenic snps located in the pericentromere regions corresponding to the ref'''
for file in $div_files
do
#query=Oh43
query=$(echo $file | cut -f1 -d "." | cut -f2 -d "_")
pericentromere_file=${ref}_${query}.pericentromere.aligned.div.txt
bedtools closest -a ${ref}.pericentromere.bed -b ${ref}_${query}.aligned.div.txt > $pericentromere_file
#start_coord=$(cat ${ref}.pericentromere.bed | cut -f2)
#end_coord=$(cat ${ref}.pericentromere.bed | cut -f3)
#syn_snp_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window_${ref}/${query}_${chr}.mapped-to.${ref}_${chr}-${ref_start_coords}_${ref_end_coords}-${ref_start_coords}_${ref_end_coords}.syn.snp.vcf
syn_snp_sum_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${query}_${chr}.mapped-to.${ref}_${chr}.pericentromere.syn.snp.bed
if [ -f "$syn_snp_sum_file" ]; then
  rm $syn_snp_sum_file
fi

while IFS= read line
do
ref_start_coords=$(echo $line | cut -f4- -d " "| awk -v var1=$ref -v var2=$query -v OFS=' ' '{if($1==var1){print$2}else{print$5}}')
ref_end_coords=$(echo $line |cut -f4- -d " "| awk -v var1=$ref -v var2=$query -v OFS=' ' '{if($1==var1){print$3}else{print$6}}')
query_start_coords=$(echo $line |cut -f4- -d " "| awk -v var1=$ref -v var2=$query -v OFS=' ' '{if($4==var2){print$5}else{print$2}}')
query_end_coords=$(echo $line |cut -f4- -d " "| awk -v var1=$ref -v var2=$query -v OFS=' ' '{if($4==var2){print$6}else{print$3}}')
echo -e $ref"\t"$ref_start_coords"\t"$ref_end_coords"\t"$query"\t"$query_start_coords"\t"$query_end_coords
syn_snp_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window_${ref}/${query}_${chr}.mapped-to.${ref}_${chr}-${ref_start_coords}_${ref_end_coords}-${ref_start_coords}_${ref_end_coords}.syn.snp.vcf
if [ -f "$syn_snp_file" ]; then
cat $syn_snp_file | \
awk -v OFS='\t' -v var1=${ref} -v var2=${query} '{print var1,var2,$0}' | \
awk -v OFS='\t' -v var1=${start_coord} -v var2=${end_coord} '{if($4>var1&&$4<var2){print$0}}' >> $syn_snp_sum_file
else
  echo "snp file not found"
fi
done < $pericentromere_file
done

'''merge snps'''
cat /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/*_${chr}.mapped-to.${ref}_${chr}.pericentromere.syn.snp.bed > ${ref}.NAM.pericentromere.syn.snp.bed

'''extract intergenic snps'''
intergenic_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency/${ref}.${chr}.nonUMR_gene.bed
cat ${ref}.NAM.pericentromere.syn.snp.bed | \
awk -v OFS='\t' '{print $3,$4,$4+1,$0}' | \
bedtools intersect -a - -b $intergenic_file -wa| \
cut -f4- > ${ref}.NAM.pericentromere.syn.intergenic.snp.bed

ref_genome=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.fasta
cat B73.NAM.pericentromere.syn.intergenic.snp.bed | cut -f4 |sort | uniq | \
awk -v OFS='\t' -v var1=$chr '{print var1,$0,$0+1}' | \
bedtools getfasta -fi $ref_genome -bed - -bedOut > ${ref}.NAM.pericentromere.syn.intergenic.allsnps.bed

lines=$(cat ${ref}.NAM.pericentromere.syn.intergenic.snp.bed |cut -f2 |sort |uniq)
cat ${ref}.NAM.pericentromere.syn.intergenic.allsnps.bed | \
awk -v OFS='\t' -v var1=$ref '{print var1,$0,$4}' > ${ref}.NAM.pericentromere.syn.intergenic.allsnps.sum
for line in $lines
do
syn_aligned_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${ref}_${line}.syn.aligned.bed
cat $syn_aligned_file | \
awk -v OFS='\t' -v var1=$chr '{print var1,$2,$3}' > ${ref}_${line}.syn.coord.tmp.bed
cat ${ref}.NAM.pericentromere.syn.intergenic.snp.bed | \
awk -v OFS='\t' -v var1=$line '{if($2==var1){print$3,$4,$4+1,$6}}' | \
bedtools intersect -b - -a ${ref}.NAM.pericentromere.syn.intergenic.allsnps.bed -wao| \
bedtools intersect -a - -b ${ref}_${line}.syn.coord.tmp.bed -wao| \
awk -v OFS='\t' -v var1=$line \
'{if($13==0){print var1,$1,$2,$3,$4,"NA"}else if($13!=0&&$9==1){print var1,$1,$2,$3,$4,$8}else if($13!=0&&$9==0){print var1,$1,$2,$3,$4,$4}}' >> ${ref}.NAM.pericentromere.syn.intergenic.allsnps.sum
rm ${ref}_${line}.syn.coord.tmp.bed
done

