wdir=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte
ref=B73
NAMline=(B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
window=(1000 5000 10000 15000 20000 50000 100000)


ml BEDTools/2.29.2-GCC-8.3.0
#identify all syntenic snps
for chr in chr{1..10}
do
  echo $chr
  out_wdir=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/$chr
  mkdir -p $out_wdir
for line in ${NAMline[@]}
do
  echo $line
  cat /scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window/${line}_${chr}.mapped-to.${ref}_${chr}-*_*-*_*.syn.snp.vcf | \
  grep -v "#" | \
  sed 's|=|\t|g' | sed 's|;|\t|g' |
  awk -v var1=$line -v var2=$chr -v OFS='\t' -v var3=$ref '{print var3,$2,$2+1,var1,$11,$11+1,var2,$4,$5}' | \
  bedtools sort -i - |uniq > ${out_wdir}/${line}.${chr}.syn.snp.bed
done
done

#map snps in a window
#test different window size (1Kb -20Kb)
for window_size in ${window[@]}
do
cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.fasta.fai | \
grep "chr" | cut -f1,2 | \
bedtools makewindows -g - -w $window_size > /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.$window_size.windows
done

cd $wdir
#map snps to different window size (1Kb -20Kb)
for window_size in ${window[@]}
do
  snp_output=NAM.syn.snp.${window_size}.bed
  if [ -f "$snp_output" ]; then
      echo "$snp_output exists."
  else
  echo $window_size
  for chr in chr{1..10}
  do
  out_wdir=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/$chr
  for line in ${NAMline[@]}
  do
  cat ${out_wdir}/${line}.${chr}.syn.snp.bed | \
  awk -v var1=$chr -v OFS='\t' '{print var1,$2,$3}' | \
  bedtools intersect -a /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.${window_size}.windows -b - -c | \
  awk -v var1=$chr -v var2=$line  -v OFS='\t' '{if($1==var1){print var2, $0}}' >> NAM.syn.snp.${window_size}.bed
done
done
fi
done


#map aligned length to different window size (1Kb -20Kb)
cd $wdir
for window_size in ${window[@]}
do
  snp_aligned_output=NAM.syn.snp.alignedlen.${window_size}.bed
  if [ -f "$snp_aligned_output" ]; then
      echo "$snp_aligned_output exists."
  else
  echo $window_size
  for chr in chr{1..10}
  do
    echo $chr
  for line in ${NAMline[@]}
  do
    #echo $line
  cat NAM.syn.snp.${window_size}.bed | \
  awk -v var1=$line -v var2=$chr -v var3=$ref -v OFS='\t' '{if($1==var1&&$2==var2){print var3,$3,$4,$5}}' | \
  bedtools intersect -a - -b /scratch/jl03308/NAM_pancentromere/NAM_SV/${chr}/B73_${line}.aligned.bed -wao | \
  bedtools sort -i - | \
  bedtools groupby -i - -g 1,2,3,4 -c 17 -o sum | \
  awk -v var1=$line -v var2=$chr -v OFS='\t' '{if($5>0){print var1,var2,$2,$3,$4,$5,$4/$5/2/3.3 *100000000}else{print var1,var2,$2,$3,$4,$5,"NA"}}' >> NAM.syn.snp.alignedlen.${window_size}.bed
done
done
fi
ml R
Rscript ~/git/NAM_pancentromere/NAM_SV/snp_distribution.window.R $wdir NAM.syn.snp.alignedlen.${window_size}.bed NAM.syn.snp.alignedlen.${window_size}.png
done
