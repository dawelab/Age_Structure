#!/bin/bash

ref='B73'

ml R
NAMline=(B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)

cd /scratch/jl03308/NAM_pancentromere/methylation/UMRs

# if [ -f "${ref}_NAM.aligned.bed" ]; then
#   echo "NAM aligned segments created"
# else

cat meth_B73.ref_B73.UMR.bed | grep -v "scaf" | awk '{print "chr"$1"\t"$2"\t"$3"\tB73"}' > NAM.ref_B73.UMR.syn.bed

for i in "${NAMline[@]}"
do
   echo "$i"
   for chr in {1..10}
   do
   cat meth_${i}.ref_B73.UMR.bed | \
   awk -v var1=$chr -v var2=$i '{if($1==var1){print"B73\t"$2"\t"$3}}' | \
   bedtools intersect -a - -b /scratch/jl03308/NAM_pancentromere/NAM_SV/chr${chr}/B73_${i}.aligned.syn.bed | \
   awk -v var1=$chr -v var2=$i '{print "chr"var1"\t"$2"\t"$3"\t"var2}' >> NAM.ref_B73.UMR.syn.bed
   # or do whatever with individual element of the array
done
done

#fi
for chr in chr{1..10}
do
Rscript ~/git/NAM_pancentromere/UMR/umrsegments.R \
/scratch/jl03308/NAM_pancentromere/methylation/UMRs/ \
/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${ref}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff \
/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum \
/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.chrom.sizes \
/scratch/jl03308/NAM_pancentromere/methylation/UMRs/NAM.ref_B73.UMR.syn.bed \
${chr} ${ref} NAM
done

chr=chr2
start=230000000
end=240000000

Rscript ~/git/NAM_pancentromere/UMR/umrsegments_snpdating.zoomin.R \
/scratch/jl03308/NAM_pancentromere/methylation/UMRs/ \
/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${ref}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff \
/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum \
/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.chrom.sizes \
/scratch/jl03308/NAM_pancentromere/methylation/UMRs/NAM.ref_B73.UMR.syn.bed \
${chr} ${ref} NAM ${start} ${end}

# umrs mappable to B73 but variable across NAM
cat NAM.ref_B73.UMR.syn.bed | bedtools sort -i - | bedtools merge -i - >  NAM.ref_B73.mergedUMR.syn.bed

ml BEDTools/2.29.2-GCC-8.3.0
cd /scratch/jl03308/NAM_pancentromere/methylation/UMRs
for i in "${NAMline[@]}"
do
  echo "$i"
  for chr in {1..10}
  do
   cat NAM.ref_B73.UMR.syn.bed | \
   awk -v var1=$chr -v var2=$i '{if($1==var1&&$4==var2){print"B73\t"$2"\t"$3}}' | \
   awk -v var1=$chr '{if($1==var1){print var1"\t"$2"\t"$3}}' | \
   bedtools intersect -a - -b NAM.ref_B73.mergedUMR.bed -v | \
   awk -v var1=$i '{print var1"\t"$2"\t"$3}' >> NAM.ref_B73.UMR.geneticbackground.bed
 done
done
