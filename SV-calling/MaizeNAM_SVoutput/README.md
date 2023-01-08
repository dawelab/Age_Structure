## Pairwise alignment among NAM lines across 10 chromosomes  
#### Files could be used for similarity/dissimilarity calculation, classfication of NAM lines in defined regions, and TE-SV comparisons.
### syn.bed
- Pairwise syntenically aligned regions across NAM lines (Output of alignment chaining; 325 pairwise comparisons for each chromosome) are in files: NAM.chrx.pairwise.syn.bed.gz 
### sv.bed
- Tandem duplications, inversions and translocations (Output of SV characterization; 325 pairwise comparisons for each chromosome) are in files: NAM.chrx.pairwise.sv.bed.gz 

#### The cols in NAM.chrx.pairwise.syn.bed.gz are: ref_genome, ref_start,ref_end, query_genome, query_start,query_end, alignment_orientation, alignment_length. The cols in NAM.chrx.pairwise.sv.bed.gz are: ref_genome, ref_start,ref_end, query_genome, query_start,query_end, sv_type (tandemdup_q indicates tandem duplication in query while tandemdup_r represents tandem duplication in reference genome).

