Classify UMRs into syntenic UMRs and neoUMRs based on their genetic background compared with B73
```
NAMline=(B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
div_file=/scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergeni_div.bed
cd /scratch/jl03308/NAM_pancentromere/methylation/UMRs/UMRs_classification
for line in "${NAMline[@]}"
do
cat $div_file | \
awk -v var2=$line -v OFS="\t" '{if($4==var2){print $10,$5,$6,$9}}' | \
bedtools sort -i - | \
sed 's|chr||g' | \
bedtools intersect -b - -a /scratch/jl03308/NAM_pancentromere/methylation/UMRs/meth_${line}.ref_${line}.UMR.bed -wao -f 0.5| \
awk -v OFS='\t' '{if($11>0){print $1,$2,$3,$5,$10}else{print $1,$2,$3,$5,"unaligned"}}' | \
bedtools sort -i - | \
grep -v "scaf" | \
awk -v OFS='\t' -v var1=$chr -v var2=$line '{print var2, "chr"$0}'  >> NAM.line_UMR.classification.bed
done

```

Identify the genes that individual UMRs regulate across NAM lines based on the condition that -10bp to +400bp from gene transcription start site overlap with UMR
```
for line in "${NAMline[@]}"
do
cat NAM.line_UMR.classification.bed | \
awk -v var2=$line -v OFS="\t" '{if($1==var2){print $2,$3,$4,$5,$6}}' | \
bedtools sort -i - > $line.line_UMR.classification.bed

cat $line.TSSrange.bed | \
awk -v var1=$chr -v var2=$line -v OFS="\t" '{print $1,$4,$5,$7,$9}' | \
grep -v "scaf" | \
bedtools sort -i - | \
bedtools intersect -b - -a $line.line_UMR.classification.bed -wao | \
awk -v var2=$line -v OFS="\t" '{if($11>0){print var2,$1,$2,$3,$4,$5,$10,$9}}' >> NAM.line_UMR.TSSrange.bed
done
```

Calculate ifferential expression value for each syntenic gene pair between NAM and B73

`umr_diffexpression.R` 
