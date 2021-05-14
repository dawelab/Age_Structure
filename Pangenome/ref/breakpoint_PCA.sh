'''B73 as reference'''
ref=B73
NAMline=(B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
#echo $(shuf -e "${NAMline[@]}")
'''merge alignment files'''
for chr in chr{1..10}
do
dir=/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}
for query in "${NAMline[@]}"
do
  if [ "$query" != "$ref" ]
  then
  cat $dir/${ref}_${query}.syn.aligned.bed | \
  cut -f1-10 | \
  awk -v var1=$chr '{print $0"\t"var1}' >> /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${ref}.syn.aligned.bed
fi
done
done

'''count alignment frequency, alignment permutation'''
ml Anaconda3
python ~/git/NAM_pancentromere/pangenome_ref/allele_frequency_cal.py \
-i /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${ref}.syn.aligned.bed \
-g /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/NAM.genome \
-r $ref \
-d /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis \
-f /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/alignment_frequency_allchrs.bed \
-s /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/alignment_permutation.bed


'''gene/UMR analysis based on alignments'''
for chr in chr{1..10}
do
wdir=/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency
mkdir -p $wdir
cat /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/alignment_frequency_allchrs.bed | \
awk -v var1=$chr -v OFS='\t' '{if($28==var1){print$1,$2,$29}}' > $wdir/$chr.allele_frequency.bed
sh ~/git/NAM_pancentromere/pangenome_ref/Gene_UMR_cal.sh $chr $ref
done
cat $wdir/${ref}.chr*.stat > /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/$ref.allele_frequency.stat

'''plot alignment frequency (bar plot)'''
ml R
Rscript /home/jl03308/git/NAM_pancentromere/pangenome_ref/plot_allele_frequency.R

'''PCA analysis of breakpoints'''
ml Anaconda3
python ~/git/NAM_pancentromere/pangenome_ref/PCA.py \
-i /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${ref}.syn.aligned.bed \
-d /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis


'''saturation analysis'''
for chr in chr{1..10}
do
  cat $chr.alignment_saturation.bed | awk -v var1=$chr '{print var1"\t"$0}' >> B73.alignment_saturation.bed
done
'''plot alignment frequency (bar plot)'''
ml R
Rscript /home/jl03308/git/NAM_pancentromere/pangenome_ref/plot_saturation.R
