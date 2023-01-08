## Table of Contents
1. [Alignment](#alignment)
2. [Chaining](#chaning)
3. [SV characterization](#sv-characterization)
4. [Plot](#plot)
5. [File format](#file-format)

## Tools used for SV calling and alignment plotting
- Alignment
    * minimap2 (2.17)
- Chaining and SV calling
    * Anaconda3 (2020.02)
- Plot
    * karyoploteR
    * ggplot2

## General procedure for pairwise genome comparison

## Alignment
Step1: minimap2 pairwise alignment:
ref, query are input fa sequences, and chr is chromosome number, output is $paf_file
```
sh minimap2_alignment.sh $ref $query $chr
```
Step2: remove alignment sequences from minimap2 output file:
ref, query are input genome names
```
cat $paf_file | \
awk -v OFS='\t' -v var1=$ref -v var2=$query '{print var2"_"$1,$2,$3,$4,$5,var1"_"$6,$7,$8,$9,$10,$11,$12}' > $paf_noseq
```
## Chaining
output (aligned_file) is syntenic aligned segments after chaining, which also include alignable inversions and translocations
```
python chaining.py \
 -i $paf_noseq \
 -o $aligned_file
```
## SV characterization
aligned_file is output of chaining, selected_aligned_file is upon filtering inversions/tranlocations by a 20Kb cutoff,
ref_chrsize and query_chrsize are chromosome size (bp), sv_file is the output file which summarizes sv and rearrangements,
unaligned_file is the pairwise unalignment result (SV final result), fig_file summarizes the SV distribution
```
python sv_detect.py \
-i $aligned_file -a $selected_aligned_file -gs1 $ref_chrsize -gs2 $query_chrsize -s $sv_file -o $unaligned_file -f $fig_file
```
## Plot
Pairwise alignement plot:
aligned_file is output of chaining, centromere_file is a bed file that marks the location of centromeres,
gff_ref and gff_query are TE annotation files, gap_file is a bed file that marks the location of gaps,
ref, query are input genome names, and chr is chromosome number
```
Rscript pairwise_alignedsegments.karyoploter.R \
 $out_dir \
 $paf_noseq \
 $aligned_file \
 $centromere_file \
 $gff_ref \
 $gff_query \
 $gap_file \
 $ref $query $chr
```
## File format
Pairwise syntenically aligned regions across NAM lines (Output of alignment chaining; 325 pairwise comparisons for each chromosome) are in files: NAM.chrx.pairwise.syn.bed.gz 
Tandem duplications, inversions and translocations (Output of SV characterization; 325 pairwise comparisons for each chromosome) are in files: NAM.chrx.pairwise.sv.bed.gz 
The cols in NAM.chrx.pairwise.syn.bed.gz are: ref_genome, ref_start,ref_end, query_genome, query_start,query_end, alignment_orientation, alignment_length
The cols in NAM.chrx.pairwise.sv.bed.gz are: ref_genome, ref_start,ref_end, query_genome, query_start,query_end, sv_type (tandemdup_q indicates tandem duplication in query while tandemdup_r represents tandem duplication in reference genome).

