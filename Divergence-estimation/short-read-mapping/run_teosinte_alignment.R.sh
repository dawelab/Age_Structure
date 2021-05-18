
ref='B73'
query='P39'
chr='chr3'

wdir=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun
gff_file1=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${ref}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff
gff_file2=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${query}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff
centromere_file=/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum
genome_size_file=${ref}.*pseudomolecules-v*.chrom.sizes
paf_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${query}_${chr}.mapped-to.${ref}_${chr}.sorted.noseq.paf
gap_file=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/NAM.100ngap.bed
alignment_file=/scratch/jl03308/NAM_pancentromere/NAM_SV/${chr}/${ref}_${query}.aligned.bed

Rscript ~/git/NAM_pancentromere/teosinte/pairwise_alignedsegments.teosinte_readalignment.R \
$wdir \
$paf_file \
$centromere_file \
$gff_file1 \
$gff_file2 \
$gap_file \
$ref \
$query \
$chr \
$alignment_file
