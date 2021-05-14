ml BEDTools/2.29.2-GCC-8.3.0
ref=B73
NAMline=(B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)

# chr=chr9
# repeat_start=22
# repeat_end=1744680


for syn_repeat_file in /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/arrays/*knob*.bed
 do
   chr=$(basename $syn_repeat_file | cut -f1 -d ".")
   prefix=$(basename $syn_repeat_file | cut -f2 -d ".")
   repeat_type=$(basename $syn_repeat_file | cut -f2 -d "." |sed 's/[0-9]*//g')
   repeat_start=$(cat $syn_repeat_file |grep $ref | cut -f2)
   repeat_end=$(cat $syn_repeat_file |grep $ref | cut -f3)
   if [[ $repeat_end ]]
   then
   wdir=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/$chr
   mkdir -p $wdir
   #extract monomer location in ref coordinates
   cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${ref}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff | \
awk '{if($3~"knob"||$3~"Cent"||$3~"rDNA"||$3~"subtelomere"){print$0}}' | \
cut -f1,4,5 |cut -f2 -d "_" | \
awk -v var1=${chr} -v OFS='\t' '{if($1==var1){print "B73",$2,$3}}' | \
bedtools merge -i - -d 10 > ${wdir}/${ref}.${chr}.${prefix}.monomer.bed

for line in ${NAMline[@]}
do
  #extract syn snps of each line in knob
cat /scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/${chr}/${line}.${chr}.syn.snp.bed | \
awk -v var1=$repeat_start -v var2=$repeat_end '{if($2>=var1&&$3<=var2){print$0}}' | \
bedtools intersect -b ${wdir}/${ref}.${chr}.${prefix}.monomer.bed -a - -v -wb > ${wdir}/${chr}.${line}.${prefix}.syn.snp.no-monomer.bed
  #extract knob regions that are TEs (no-monomer), map syn snps to TEs, calculate age
cat /scratch/jl03308/NAM_pancentromere/NAM_SV/${chr}/${ref}_${line}.aligned.syn.bed | \
awk -v var1=$repeat_start -v var2=$repeat_end '{if($2>=var1&&$3<=var2){print$0}}' | \
bedtools subtract -a - -b ${wdir}/${ref}.${chr}.${prefix}.monomer.bed | \
bedtools intersect -a - -b ${wdir}/${chr}.${line}.${prefix}.syn.snp.no-monomer.bed -c | \
awk -v OFS='\t' '{print $0, $3-$2+1, $13/($3-$2+1)/2/3.3*100000000}'  >> ${wdir}/${chr}.${prefix}.syn.snp.dating.bed
done
fi
done

cd /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/
for file in chr*/chr*.knob*.syn.snp.dating.bed
do
  prefix=$(basename $file | cut -f1,2 -d ".")
  cat $file | awk -v OFS='\t' -v var1=$prefix '{print $0,var1}' >> NAM.knob.dating.bed
done
