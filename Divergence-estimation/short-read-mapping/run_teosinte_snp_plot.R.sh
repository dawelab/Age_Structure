ref=B73

ml BEDTools/2.29.2-GCC-8.3.0
ref_index=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.fasta.fai
if [ -f "$ref_index" ]; then
    echo "$ref_index exists."
else
  ml SAMtools/1.9-GCC-8.3.0
  samtools faidx /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.fasta
fi

cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.fasta.fai | \
grep "chr" | cut -f1,2 | \
bedtools makewindows -g - -w 10000 > /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.10Kb.windows

wdir=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte
ml R
ml BEDTools/2.29.2-GCC-8.3.0
lines=($(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centC.array.assembly.status.content | cut -f1 |sort |uniq))
teolines=(TIL03 TIL04 TIL05 TIL06 TIL09 TIL11 TIL14 TIL15 Pav_ISU_1 Mex_ISU_1 Mex_ISU_2 TIL08 TIL25 Hue-2 Hue-4 dip1 Trip_ISU_1 TDD39103)

for line in "${lines[@]}"
do
  cat /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.${ref}.flt.snps.bed | \
  grep -v "scaf" | \
  bedtools intersect -a /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.10Kb.windows -b - -c | \
  awk -v var1=$line '{print var1"\t"$0"\tNAM"}' >> ${wdir}/NAM.${ref}.flt.vcf.bed
done

for line in "${teolines[@]}"
do
 cat /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.${ref}.flt.snps.bed | \
 grep -v "scaf" | \
 bedtools intersect -a /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.10Kb.windows -b - -c | \
 awk -v var1=$line '{print var1"\t"$0"\tTeo"}' >> ${wdir}/Teo.${ref}.flt.vcf.bed
done

cat ${wdir}/Teo.${ref}.flt.vcf.bed ${wdir}/NAM.${ref}.flt.vcf.bed |grep -v $ref > ${wdir}/Teo_NAM.${ref}.flt.vcf.bed

ref_gff=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${ref}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff
ref_genomesize=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.chrom.sizes
if [ -f "$ref_genomesize" ]; then
    echo "$ref_genomesize exists."
else
  cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.fasta.fai | cut -f1,2 |grep "chr" > $ref_genomesize
fi
centromere_file=/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum
output_file=${wdir}/Teo_NAM.${ref}.flt.vcf.bed

for chr in chr{1..10}
do
Rscript /home/jl03308/git/NAM_pancentromere/teosinte/snp_plot.R \
$wdir \
$ref_gff  \
$centromere_file \
$ref_genomesize \
$output_file \
${ref} \
${chr} \
Teo_NAM_${ref}
done

pc_parvlines=($(ls -d -- /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/PC* | cut -f8 -d "/"))
for line in "${pc_parvlines[@]}"
do
 cat /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.${ref}.flt.snps.bed | \
 grep -v "scaf" | \
 bedtools intersect -a /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.10Kb.windows -b - -c | \
 awk -v var1=$line '{print var1"\t"$0"\tTeo"}' >> ${wdir}/PC.${ref}.flt.vcf.bed
done

#sample N11 does not map to maize genomes
cat ${wdir}/PC.${ref}.flt.vcf.bed | grep -v "N11" | cat - ${wdir}/NAM.${ref}.flt.vcf.bed | grep -v $ref > ${wdir}/PC_NAM.${ref}.flt.vcf.bed
output_file=${wdir}/PC_NAM.${ref}.flt.vcf.bed

for chr in chr{1..10}
do
Rscript /home/jl03308/git/NAM_pancentromere/teosinte/snp_plot.PC.R \
$wdir \
$ref_gff  \
$centromere_file \
$ref_genomesize \
$output_file \
${ref} \
${chr} \
PC_${ref}
done


wdir=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte
ml R
ml BEDTools/2.29.2-GCC-8.3.0
lines=($(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centC.array.assembly.status.content | cut -f1 |sort |uniq))
teolines=(TIL03 TIL04 TIL05 TIL06 TIL09 TIL11 TIL14 TIL15 Pav_ISU_1 Mex_ISU_1 Mex_ISU_2 TIL08 TIL25 Hue-2 Hue-4 dip1 Trip_ISU_1 TDD39103)


'''snp normalization by alignment gaps'''
for line in "${lines[@]}"
do
  cat /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.q20.RPKM.bedgraph | \
  awk '{if($4>=0.25&&$4<=4){print$0}}' | \
  bedtools intersect -a /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.10Kb.windows -b - -wao | \
  bedtools sort -i - | \
  bedtools groupby -i - -g 1,2,3 -c 8 -o sum > /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.${ref}.normalized10Kb.windows
  cat /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.q20.flt.snps.bed | \
  grep -v "scaf" | \
  bedtools intersect -a /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.${ref}.normalized10Kb.windows -b - -c | \
  awk -v var1=$line -v OFS='\t' '{if($4!=0){print var1,$1,$2,$3,$5/$4*10000,"NAM",$4}else{print var1,$1,$2,$3,$5,"NAM",$4}}' >> ${wdir}/NAM.${ref}.flt.vcf.normalized.bed
done

for line in "${teolines[@]}"
do
  cat /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.${ref}.RPKM.bedgraph| \
  awk '{if($4>=0.25&&$4<=4){print$0}}' | \
  bedtools intersect -a /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.10Kb.windows -b - -wao | \
  bedtools sort -i - | \
  bedtools groupby -i - -g 1,2,3 -c 8 -o sum > /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.${ref}.normalized10Kb.windows
  cat /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.${ref}.flt.snps.bed | \
  grep -v "scaf" | \
  bedtools intersect -a /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.${ref}.normalized10Kb.windows -b - -c | \
  awk -v var1=$line -v OFS='\t' '{if($4!=0){print var1,$1,$2,$3,$5/$4*10000,"Teo",$4}else{print var1,$1,$2,$3,$5,"Teo",$4}}' >> ${wdir}/Teo.${ref}.flt.vcf.normalized.bed
done

cat ${wdir}/Teo.${ref}.flt.vcf.normalized.bed ${wdir}/NAM.${ref}.flt.vcf.normalized.bed|grep -v $ref > ${wdir}/Teo_NAM.${ref}.flt.vcf.normalized.bed

ref_gff=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${ref}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff
ref_genomesize=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.chrom.sizes
if [ -f "$ref_genomesize" ]; then
    echo "$ref_genomesize exists."
else
  cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.fasta.fai | cut -f1,2 |grep "chr" > $ref_genomesize
fi
centromere_file=/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum
output_file=${wdir}/Teo_NAM.${ref}.flt.vcf.normalized.bed

for chr in chr{1..10}
do
Rscript /home/jl03308/git/NAM_pancentromere/teosinte/snp_plot.normalized.R \
$wdir \
$ref_gff  \
$centromere_file \
$ref_genomesize \
$output_file \
${ref} \
${chr} \
Teo_NAM_${ref}
done


# for line in "${lines[@]}"
# do
#   cat ${wdir}/NAM.${ref}.flt.vcf.bed | awk -v var1=$line -v OFS='\t' '{if($1==var1){print $2,$3,$4,$5,$6}}' | \
#   grep -v "scaf" | \
#   bedtools intersect -a /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.${ref}.normalized10Kb.windows -b - -wao | \
#   awk -v var1=$line -v OFS='\t' '{if($4!=0){print var1,$1,$2,$3,$8/$4*10000,"NAM",$4}else{print var1,$1,$2,$3,$8,"NAM",$4}}' >> ${wdir}/NAM.${ref}.flt.vcf.normalized.bed
# done
