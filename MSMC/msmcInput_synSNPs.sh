#!/bin/bash

chr="$1"

ml msmc/2.1.2-GCC-8.3.0
ml Anaconda3
ml BCFtools/1.6-foss-2019b
ml BEDTools/2.29.2-GCC-8.3.0

wdir=/scratch/jl03308/NAM_pancentromere/teosinte_divergence/MSMC/Syn/$chr
mkdir -p $wdir
nam=(B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
ref=B73
line1=B97
line2=CML103
line3=CML228
line4=CML247
line5=CML277
line6=CML322
line7=CML333
line8=CML52
line9=CML69
line10=HP301
line11=IL14H
line12=Ki11
line13=Ki3
line14=Ky21
line15=M162W
line16=M37W
line17=Mo18W
line18=MS71
line19=NC350
line20=NC358
line21=Oh43
line22=Oh7b
line23=P39
line24=Tx303
line25=Tzi8

vcf_dir=/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment

for line in "${nam[@]}"
do
#cat /scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window/${line}_${chr}.mapped-to.${ref}_${chr}-*_*-*_*.syn.snp.vcf | sort -k2,2n > ${vcf_dir}/${line}/${line}.${chr}.syn.vcf
cat ${vcf_dir}/${line}/${line}.${chr}.syn.vcf | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,"PASS", ".",$9,$10}'  > ${vcf_dir}/${line}/${line}.${chr}.syn.format.vcf
zcat ${vcf_dir}/${line}/${line}.q20.flt.vcf.gz | grep "#"  |cat - ${vcf_dir}/${line}/${line}.${chr}.syn.format.vcf > ${vcf_dir}/${line}/${line}.${chr}.syn.header.vcf
zippedfile=${vcf_dir}/${line}/${line}.${chr}.syn.header.vcf.gz
if [[ -f "$zippedfile" ]]; then
  rm $zippedfile
fi
bgzip ${vcf_dir}/${line}/${line}.${chr}.syn.header.vcf
bcftools index -f ${vcf_dir}/${line}/${line}.${chr}.syn.header.vcf.gz
# cat ${vcf_dir}/${line}/${line}.${chr}.syn.vcf | awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$7,"PASS", $9, ".", "1|1"}' > ${vcf_dir}/${line}/${line}.${chr}.syn.phased.vcf
# zcat ${vcf_dir}/${line}/${line}.${chr}.syn.header.vcf.gz | grep "#" | tail -1 |cat - ${vcf_dir}/${line}/${line}.${chr}.syn.phased.vcf > ${vcf_dir}/${line}/${line}.${chr}.syn.phased.header.vcf
# zippedfile=${vcf_dir}/${line}/${line}.${chr}.syn.phased.header.vcf.gz
# if [[ -f "$zippedfile" ]]; then
#   rm $zippedfile
# fi
# bgzip ${vcf_dir}/${line}/${line}.${chr}.syn.phased.header.vcf
done

# # mask bed
# for line in "${nam[@]}"
# do
#   cat /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/B73_${line}.syn.aligned.bed | \
#   awk -v var1=$chr '{print var1"\t"$2"\t"$3}' | \
#   bedtools merge -i - > /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line}.mask.bed
# done

bcftools merge ${vcf_dir}/${line1}/${line1}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line2}/${line2}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line3}/${line3}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line4}/${line4}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line5}/${line5}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line6}/${line6}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line7}/${line7}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line8}/${line8}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line9}/${line9}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line10}/${line10}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line11}/${line11}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line12}/${line12}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line13}/${line13}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line14}/${line14}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line15}/${line15}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line16}/${line16}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line17}/${line17}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line18}/${line18}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line19}/${line19}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line20}/${line20}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line21}/${line21}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line22}/${line22}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line23}/${line23}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line24}/${line24}.${chr}.syn.header.vcf.gz \
${vcf_dir}/${line25}/${line25}.${chr}.syn.header.vcf.gz \
-Oz --force-samples > ${wdir}/NAM.${chr}.syn.header.vcf.gz

#remove duplicates
bcftools norm ${wdir}/NAM.${chr}.syn.header.vcf.gz -d both -Oz -o ${wdir}/NAM.${chr}.syn.header.norm.vcf.gz
#phase snp calls across NAM lines
ml Java/1.8.0_241
java -Xmx240g -jar ~/beagle.18May20.d20.jar gt=${wdir}/NAM.${chr}.syn.header.norm.vcf.gz chrom=$chr impute=false out=${wdir}/NAM.${chr}.syn.phased

# extract phased haplotypes of each line
ml VCFtools/0.1.16-GCC-8.3.0-Perl-5.30.0
for line in "${nam[@]}"
do
vcftools --gzvcf ${wdir}/NAM.${chr}.syn.phased.vcf.gz --indv /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${line}/${line}.q20.sorted.bam --recode --out ${wdir}/${line}.${chr}.syn.phased
bgzip ${wdir}/${line}.${chr}.syn.phased.recode.vcf
done

ml msmc/2.1.2-GCC-8.3.0
ml Anaconda3
python ~/msmc-tools/generate_multihetsep.py --chr $chr \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line1}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line2}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line3}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line4}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line5}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line6}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line7}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line8}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line9}.mask.bed  \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line10}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line11}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line12}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line13}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line14}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line15}.mask.bed  \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line16}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line17}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line18}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line19}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line20}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line21}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line22}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line23}.mask.bed  \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line24}.mask.bed \
--mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line25}.mask.bed \
${wdir}/${line1}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line2}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line3}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line4}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line5}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line6}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line7}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line8}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line9}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line10}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line11}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line12}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line13}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line14}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line15}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line16}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line17}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line18}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line19}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line20}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line21}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line22}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line23}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line24}.${chr}.syn.phased.recode.vcf.gz \
${wdir}/${line25}.${chr}.syn.phased.recode.vcf.gz > ${wdir}/NAM.$chr.multihetsep.txt

cat ${wdir}/NAM.$chr.multihetsep.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("","\t",$4)}1' | \
cut -f1,2,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53 | \
awk '{print $1"\t"$2"\t"$3"\t"$4""$5""$6""$7""$8""$9""$10""$11""$12""$13""$14""$15""$16""$17""$18""$19""$20""$21""$22""$23""$24""$25""$26""$27""$28}'> ${wdir}/NAM.$chr.multihetsep.half.txt
msmc2 -m 0.0012 -p 5*4+25*2+5*4 -o ${wdir}/NAM.${chr}.msmc2 ${wdir}/NAM.$chr.multihetsep.half.txt


# python ~/msmc-tools/generate_multihetsep.py --chr $chr \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line1}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line2}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line3}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line4}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line5}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line6}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line7}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line8}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line9}.mask.bed  \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line10}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line11}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line12}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line13}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line14}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line15}.mask.bed  \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line16}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line17}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line18}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line19}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line20}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line21}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line22}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line23}.mask.bed  \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line24}.mask.bed \
#  --mask /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${line25}.mask.bed \
#     ${vcf_dir}/${line1}/${line1}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line2}/${line2}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line3}/${line3}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line4}/${line4}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line5}/${line5}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line6}/${line6}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line7}/${line7}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line8}/${line8}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line9}/${line9}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line10}/${line10}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line11}/${line11}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line12}/${line12}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line13}/${line13}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line14}/${line14}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line15}/${line15}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line16}/${line16}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line17}/${line17}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line18}/${line18}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line19}/${line19}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line20}/${line20}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line21}/${line21}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line22}/${line22}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line23}/${line23}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line24}/${line24}.${chr}.syn.phased.header.vcf.gz \
#     ${vcf_dir}/${line25}/${line25}.${chr}.syn.phased.header.vcf.gz > ${wdir}/NAM.$chr.multihetsep.txt
#
# cat ${wdir}/NAM.$chr.multihetsep.txt | awk 'BEGIN{FS=OFS="\t"}{gsub("","\t",$4)}1' | \
# cut -f1,2,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53 | \
# awk '{print $1"\t"$2"\t"$3"\t"$4""$5""$6""$7""$8""$9""$10""$11""$12""$13""$14""$15""$16""$17""$18""$19""$20""$21""$22""$23""$24""$25""$26""$27""$28}'> ${wdir}/NAM.$chr.multihetsep.half.txt
#
# msmc2 -m 0.0012 -p 5*4+25*2+5*4 -o ${wdir}/NAM.$chr.msmc2 ${wdir}/NAM.$chr.multihetsep.half.txt
#
