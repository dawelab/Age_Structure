## Tools used for SV calling and alignment plotting
- Pangenome calculation
    * Anaconda3 (2020.02)
- Plot
    * ggplot2

## Step1: combine and label chromosome size of all 26 NAM lines in a single file:
```
NAMline=(B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
for chr in chr{1..10}
do
for line in "${NAMline[@]}"
do
size=$(cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${line}.*pseudomolecules-v*.fasta.fai | awk -v var1=$chr '{if($1==var1){print$2}}')
echo -e $chr"\t"$line"\t1\t"$size >> /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/NAM.genome
done
done
```
## Step2: pangenome-permutation analysis:
Order of NAM lines was shuffled for a 1000 times, and pan-genome was calculated for every case.
run_pangenome_cal_shuf.sh executes `pangenome_cal_shuf_parallel.py`, which calculates pangenome space.
```
for chr in chr{1..10}
do
  sh run_pangenome_cal_shuf.sh $chr
done
```
## Step3: combine pangenome data and plot:
```
for chr in chr{1..10}
do
  cat $chr.pangenome.bed | awk -v var1=$chr '{print $0"\t"var1}' >> wholegenome.pangenome.bed
done

ml R
Rscript ~/git/NAM_pancentromere/pangenome/plot_pangenome_shuf.R
```
