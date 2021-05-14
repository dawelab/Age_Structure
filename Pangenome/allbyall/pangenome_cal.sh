ml BEDTools/2.29.2-GCC-8.3.0

NAMline=(B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
for chr in chr{1..10}
do
for line in "${NAMline[@]}"
do
size=$(cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${line}.*pseudomolecules-v*.fasta.fai | awk -v var1=$chr '{if($1==var1){print$2}}')
echo -e $chr"\t"$line"\t1\t"$size >> /lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/NAM.genome
done
done

'''pangenome-permutation analysis'''
for chr in chr{1..10}
do
  sh /home/jl03308/git/NAM_pancentromere/pangenome/run_pangenome_cal_shuf.sh $chr
done

'''plot pan-genome permutation result'''

for chr in chr{1..10}
do
  cat $chr.pangenome.bed | awk -v var1=$chr '{print $0"\t"var1}' >> wholegenome.pangenome.bed
done

ml R
Rscript ~/git/NAM_pancentromere/pangenome/plot_pangenome_shuf.R


  
