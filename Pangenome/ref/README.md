## Tools used for Pangenome calculation and plotting
- Pangenome calculation
    * Anaconda3 (2020.02)
- Plot
    * ggplot2

## Frequency of B73 genome space in NAM population
Step1: merge alignment files
```
ref=B73
NAMline=(B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
#echo $(shuf -e "${NAMline[@]}")
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
```
Step2: count alignment frequency, alignment permutation
```
ref=B73
ml Anaconda3
python allele_frequency_cal.py \
-i /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${ref}.syn.aligned.bed \
-g /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/NAM.genome \
-r $ref \
-d /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis \
-f /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/alignment_frequency_allchrs.bed \
-s /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/alignment_permutation.bed
```

## Categorize each segment by Gene/UMR overlaps

```
for chr in chr{1..10}
do
wdir=/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency
mkdir -p $wdir
cat /scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/alignment_frequency_allchrs.bed | \
awk -v var1=$chr -v OFS='\t' '{if($28==var1){print$1,$2,$29}}' > $wdir/$chr.allele_frequency.bed
sh Gene_UMR_cal.sh $chr $ref
done
cat $wdir/${ref}.chr*.stat > /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/$ref.allele_frequency.stat
```

## Plot

```
ml R
Rscript plot_allele_frequency.R
```
