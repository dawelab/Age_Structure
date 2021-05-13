#!/bin/bash

chr="$1"
ref="$2"
query="$3"

ml Anaconda3/2020.02

# extract align coord info from paf file
paf_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${query}_${chr}*.mapped-to.${ref}_${chr}.paf
if [ ! -f $paf_file ]; then
  paf_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${query}_${chr}*.mapped-to.${ref}_${chr}.sorted.paf
fi
prefix=$(basename $paf_file | cut -f1,2,3 -d ".")

out_dir=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}
mkdir -p $out_dir

cat $paf_file | \
awk -v OFS='\t' -v var1=$ref -v var2=$query '{print var2"_"$1,$2,$3,$4,$5,var1"_"$6,$7,$8,$9,$10,$11,$12}' > /scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${prefix}.noseq.paf

paf=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${prefix}.noseq.paf

# assign output files
aligned_file=${out_dir}/${ref}_${query}.aligned.bed
selected_aligned_file=${out_dir}/${ref}_${query}.syn.aligned.bed

unaligned_file=${out_dir}/${ref}_${query}.unaligned.bed
sv_file=${out_dir}/${ref}_${query}.sv.sum
fig_file=${out_dir}/${ref}_${query}.stats.pdf
ref_chrsize=$(cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.fasta.fai | awk -v var1=$chr '{if($1==var1){print$2}}')
query_chrsize=$(cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${query}.*pseudomolecules-v*.fasta.fai | awk -v var1=$chr '{if($1==var1){print$2}}')

if test -f "$aligned_file"; then
   echo "$aligned_file exists."
else

# chaining
python /home/jl03308/git/NAM_pancentromere/NAM_SV/chaining.py \
 -i $paf \
 -o $aligned_file

 # SV detect
python /home/jl03308/git/NAM_pancentromere/NAM_SV/sv_detect.py \
-i $aligned_file -a $selected_aligned_file -gs1 $ref_chrsize -gs2 $query_chrsize -s $sv_file -o $unaligned_file -f $fig_file

# plot pairwise alignment

ml R

gff_ref=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${ref}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff
gff_query=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${query}.pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff
centromere_file=/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum
gap_file=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/NAM.100ngap.bed

Rscript /home/jl03308/git/NAM_pancentromere/NAM_SV/pairwise_alignedsegments.karyoploter.R \
 $out_dir \
 $paf \
 $aligned_file \
 $centromere_file \
 $gff_ref \
 $gff_query \
 $gap_file \
 $ref $query $chr
fi
