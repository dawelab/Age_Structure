cd /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/arrays
ml BEDTools/2.29.2-GCC-8.3.0
# lines=$(cat 100Kb_knob_B73_synteny_chr.bed |cut -f1 |sort |uniq)
# for line in $lines
# do
#   cat 100Kb_knob_B73_synteny_chr.bed | \
#   awk -v OFS='\t' -v var1=$line '{if($1==var1){print $2,$3,$4,$5,$6,$8}}' | \
#   bedtools sort -i - | \
#   bedtools merge -i - -o collapse -c 6 -d 5000000 | \
#   awk -v var1=$line '{print var1"\t"$0}' >> 100Kb_knob_merged.bed
# done
#
# cat 100Kb_knob_merged.bed |sort -k2,2 -k3,3n > 100Kb_knob_merged.sorted.bed

#identify syntenic repeat position by checking the repeat plot for 26 genomes
cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | grep "chr9"  | awk '{if($5<10000000&&$6>50000){print$0}}'| cut -f1,2,4,5 | cut -f1 -d "." | cut -f1,3,4 |bedtools sort -i - | bedtools merge -i - -d 1000000
cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | grep "chr8" | awk '{if($4>150000000&&$5<175000000){print$0}}'| cut -f1,2,4,5 | cut -f1 -d "." | cut -f1,3,4 |bedtools sort -i - | bedtools merge -i - -d 13000000|grep -v "AB"
cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | grep "chr7" | awk '{if($4>0&&$5<10000000){print$0}}'| cut -f1,2,4,5 | cut -f1 -d "." | cut -f1,3,4 |bedtools sort -i - | bedtools merge -i - -d 1000000|grep -v "AB"
cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | grep "chr6"  | awk '{if($4>0&&$5<10000000){print$0}}'| cut -f1,2,4,5 | cut -f1 -d "." | cut -f1,3,4 |bedtools sort -i - | bedtools merge -i - -d 1000000|grep -v "AB"
cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | grep "chr6"  | awk '{if($4>10000000&&$5<30000000){print$0}}'| cut -f1,2,4,5 | cut -f1 -d "." | cut -f1,3,4 |bedtools sort -i - | bedtools merge -i - -d 10000000|grep -v "AB"|wc -l
cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | grep "chr5"  | awk '{if($4>190000000&&$5<210000000&&$6>50000){print$0}}'| cut -f1,2,4,5 | cut -f1 -d "." | cut -f1,3,4 |bedtools sort -i - | bedtools merge -i - -d 10000000|grep -v "AB"
cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | grep "chr4"  | awk '{if($4>210000000&&$5<240000000&&$6>50000){print$0}}'| cut -f1,2,4,5 | cut -f1 -d "." | cut -f1,3,4 |bedtools sort -i - | bedtools merge -i - -d 10000000|grep -v "AB"
cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | grep "chr3"  | awk '{if($4>180000000&&$5<200000000&&$6>50000){print$0}}'| cut -f1,2,4,5 | cut -f1 -d "." | cut -f1,3,4 |bedtools sort -i - | bedtools merge -i - -d 10000000|grep -v "AB"
cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | grep "chr2"  | awk '{if($4>200000000&&$5<220000000&&$6>50000){print$0}}'| cut -f1,2,4,5 | cut -f1 -d "." | cut -f1,3,4 |bedtools sort -i - | bedtools merge -i - -d 10000000|grep -v "AB"

cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | awk '{if($1=="B73"&&$6>100000){print$0}}' |grep -v "CentC"|grep "chr" |cut -f1,2,4,5 |cut -f1 -d "." > B73.100Krepeats.bed
#nor
cat All_Repeat_content_sum.csv| sed 's|,|\t|g' | grep "AF013103.1"|grep "chr6" |awk '{if($4>10000000&&$5<30000000&&$6>10000){print$0}}' |cut -f1,2,4,5 | cut -f1 -d "." | cut -f1,3,4 |bedtools sort -i - | bedtools merge -i - -d 20000000|grep -v "AB" > chr6.rDNA1.bed
#centc
repeat="Cent/CentC"
lines=$(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centC.array.assembly.status.content | cut -f1 |sort |uniq)
for chr in chr{1..10}
do
  for line in $lines
  do
    centro_start=$(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum | awk -v var1=$chr -v OFS='\t' -v var2=$line '{if($2==var1&&$1==var2){print$3}}')
    centro_end=$(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum | awk -v var1=$chr -v OFS='\t' -v var2=$line '{if($2==var1&&$1==var2){print$4}}')
    centro_lower=$(echo $centro_start - 5000000 |bc)
    centro_upper=$(echo $centro_end + 5000000 |bc)
    cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${line}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff | \
    awk -v var1=$repeat -v var3=${line}_${chr} -v OFS='\t' '{if($3==var1&&$1==var3){print$1,$4,$5}}' | \
    bedtools sort -i - | \
    bedtools merge -i - -d 1000000 | \
    awk -v var1=${line}_${chr} -v var2=$centro_lower -v var3=$centro_upper -v var4=$line -v OFS='\t' '{if($1==var1&&$2>=var2&&$3<=var3){print var4, $2,$3}}' | bedtools merge -i - -d 5000000 >> $chr.centc1.bed
    #echo -e $line"\t"$centro_lower"\t"$centro_upper >> $chr.centc1.bed
    #cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centC.array.assembly.status.content| awk -v var1=$chr -v var2=$centro_lower -v var3=$centro_upper -v var4=$line -v OFS='\t' '{if($2==var1&&$3>=var2&&$3<=var3&&$1==var4){print$0}}' |cut -f1,3,4 | bedtools sort -i - |bedtools merge -i - -d 10000000 >> $chr.centc1.bed
done
done

chr=chr1
start=184000000
end=186000000

chr=chr6
start=0
end=6000000
chr=chr6
start=14000000
end=26000000

chr=chr7
start=156000000
end=160000000

chr=chr8
start=159000000
end=163000000
chr=chr8
start=26000000
end=32000000

chr=chr9
start=0
end=2000000

chr=chr5
start=198500000
end=203500000

chr=chr6
prefix=knob2
ref=B73

wdir=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/TE_alignment/${chr}
mkdir -p $wdir
ref_gff=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${ref}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff
ref_genomesize=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.chrom.sizes
centromere_file=/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum
output_file=/scratch/jl03308/NAM_pancentromere/NAM_SV/${chr}/${ref}_NAM.div.txt
gap_file=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/NAM.100ngap.bed

repeat_array_file=/lustre2/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/arrays/${chr}.${prefix}.bed
start_range=$(cat $repeat_array_file | grep $ref |cut -f2 |awk '{print $1 - 1000000}')
end=$(cat $repeat_array_file | grep $ref |cut -f3 |awk '{print $1 + 1000000}')
if [ "$start_range" -gt 0 ]; then
  start=$start_range
else
  start=0
fi

ml R
Rscript /home/jl03308/git/NAM_pancentromere/NAM_SV/snp_dating_plot.zoomin.R \
 $wdir \
 $ref_gff  \
 $centromere_file \
 $ref_genomesize \
 $output_file \
 ${chr} ${start} ${end}

output_file=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/Teo_NAM.${ref}.flt.vcf.normalized.bed
Rscript /home/jl03308/git/NAM_pancentromere/teosinte/snp_plot.normalized.zoomin.R \
$wdir \
$ref_gff  \
$centromere_file \
$ref_genomesize \
$output_file \
${ref} \
${chr} \
Teo_NAM_${ref} ${start} ${end}


NAM_lines=$(cat $centromere_file | cut -f1 |sort |uniq |grep -v $ref)
for query in $NAM_lines
do
paf_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${query}_${chr}.mapped-to.${ref}_${chr}.*.noseq.paf
alignment_file=/scratch/jl03308/NAM_pancentromere/NAM_SV/${chr}/${ref}_${query}.aligned.bed
query_gff=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${query}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff
qstart_range=$(cat $repeat_array_file | grep $query |cut -f2 |awk '{print $1 - 1000000}')
qend=$(cat $repeat_array_file | grep $ref |cut -f3 |awk '{print $1 + 1000000}')
if [ "$start_range" -gt 0 ]; then
  qstart=$qstart_range
else
  qstart=0
fi
Rscript /home/jl03308/git/NAM_pancentromere/syn_repeat/TE_pairwisealignment.R \
 $wdir \
 $paf_file \
 $alignment_file \
 $centromere_file \
 $ref_gff  \
 $query_gff \
 $gap_file \
 ${ref} \
 ${query} \
 ${chr} ${start} ${end} ${qstart} ${qend}
done


# gff_file1 <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.gff'
# gff_file2 <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/B97.pseudomolecules-v1.fasta.mod.EDTA.TEanno.gff'
# centromere_file <- '/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum'
# genome_size_file <- 'B73.PLATINUM.pseudomolecules-v1.chrom.sizes'
# paf_file <- '/scratch/jl03308/NAM_pancentromere/genome_alignment/chr9_shujun/B97_chr9.mapped-to.B73_chr9.sorted.noseq.paf'
# #core_file <- 'B73_NAM.core.bed'
# chr <- 'chr9'
# alignment_file <- '/scratch/jl03308/NAM_pancentromere/NAM_SV/chr9/B73_B97.aligned.bed'
# gap_file <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/NAM.100ngap.bed'
# ref<-'B73'
# query<-'B97'
#
# query<-'CML103'
# gff_file2 <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/CML103.pseudomolecules-v1.fasta.mod.EDTA.TEanno.gff'
# alignment_file <- '/scratch/jl03308/NAM_pancentromere/NAM_SV/chr9/B73_CML103.aligned.bed'
#
# query<-'CML69'
# gff_file2 <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/CML69.pseudomolecules-v1.fasta.mod.EDTA.TEanno.gff'
# alignment_file <- '/scratch/jl03308/NAM_pancentromere/NAM_SV/chr9/B73_CML69.aligned.bed'
#
