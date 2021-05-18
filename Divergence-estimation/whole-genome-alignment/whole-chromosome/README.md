## Tools used for SNP calling and divergence dating
- SNP calling
    * minimap2 (v3.5)
- Divergence calculation for each syntenic aligned segment
    * Anaconda3 (2020.02)
- Extract SNPs mapped to intergenic spaces
    * bedtools (2.29.2)
- Plot
    * karyoploteR

## SNP calling
`paftools_sv_call_ref.sh` performs SNP and indel calling with paftools through processing alignments in PAF format (https://github.com/lh3/minimap2/blob/master/misc/README.md)
```
#!/bin/bash
chr="$1"
ref="$2"
synpaf_file="$3"

ref_genome=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.fasta
prefix=$(basename $synpaf_file |cut -f1-4 -d ".")
#chr=chr8
#genome_query=B97
#module load BEDTools/2.29.2-GCC-8.3.0
ml minimap2
cat $synpaf_file | sed 's|B73_||g' | sort -k8,8n | /home/jl03308/bin/k8 /home/jl03308/minimap2/misc/paftools.js call -L50 -q0 -l50 -f $ref_genome - > ${prefix}.vcf
```

## Divergence calculation for each syntenic aligned segment
`snp_dating.py` executes `paftools_sv_call_ref.sh` and summarizes SNP and divergence time for each aligned segment.

Run `snp_dating.py` for each syntenic alignments between each pair of genome.
```
#!/bin/bash
chr="$1"
ref="$2"
query="$3"
#chr=chr1
#ref=B73
#query=P39
if [ $ref = $query ]; then
  echo "pass, 100% match"
else
module load BEDTools/2.29.2-GCC-8.3.0
module load Anaconda3/5.0.1
mkdir -p /scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window_${ref}
#allpaf=$(ls /scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/*_${chr}.mapped-to.${ref}_${chr}.sorted.paf |grep -v "Ab10")

# for paf_file in $allpaf
# do
#query=$(basename $paf_file |cut -f1 -d "." |cut -f1 -d "_")
FILE=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${ref}_${query}.syn.aligned.bed

if [ -f "$FILE" ]; then
    aligned_syn_file=$FILE
    paf_unsorted_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${query}_${chr}.mapped-to.${ref}_${chr}.paf
    paf_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${query}_${chr}.mapped-to.${ref}_${chr}.sorted.paf
    cat $paf_unsorted_file | sort -k6,6 -k8,8n > $paf_file
    #only include snps of syntenic aligned segments and inversions bigger than 20Kb
    output_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${ref}_${query}.aligned.div.txt
    if [ -f "$output_file" ]; then
        echo "$output_file exists."
        rm $output_file
    fi
    while IFS= read line
    do
    ref_start_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($1==var1){print$2}else{print$5}}')
    ref_end_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($1==var1){print$3}else{print$6}}')
    query_start_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($4==var2){print$5}else{print$2}}')
    query_end_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($4==var2){print$6}else{print$3}}')
    echo -e $ref"\t"$ref_start_coords"\t"$ref_end_coords"\t"$query"\t"$query_start_coords"\t"$query_end_coords
    paf_syn_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window_${ref}/${query}_${chr}.mapped-to.${ref}_${chr}-${ref_start_coords}_${ref_end_coords}.syn.paf
    cat $paf_file | \
    awk -v var1=$ref_start_coords -v var2=$ref_end_coords '{if($8==var1&&$9==var2){print$0}}' > $paf_syn_file
    if [ -s $paf_syn_file ]
    then
    python /home/jl03308/git/NAM_pancentromere/NAM_SV/snp_dating4.py \
    -r $ref \
    -i1 $paf_syn_file \
    -i2 $aligned_syn_file \
    -s $ref_start_coords -e $ref_end_coords \
    -qs $query_start_coords -qe $query_end_coords \
    -o $output_file
    rm $paf_syn_file
    else
      echo "empty"
    fi
    done < $aligned_syn_file

else
    aligned_syn_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${query}_${ref}.syn.aligned.bed
    paf_unsorted_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${ref}_${chr}.mapped-to.${query}_${chr}.paf
    paf_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${ref}_${chr}.mapped-to.${query}_${chr}.sorted.paf
    cat $paf_unsorted_file | sort -k6,6 -k8,8n > $paf_file
    #only include snps of syntenic aligned segments and inversions bigger than 20Kb
    output_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${query}_${ref}.aligned.div.txt
    if [ -f "$output_file" ]; then
        echo "$output_file exists."
        rm $output_file
    fi
    while IFS= read line
    do
    ref_start_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($1==var1){print$2}else{print$5}}')
    ref_end_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($1==var1){print$3}else{print$6}}')
    query_start_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($4==var2){print$5}else{print$2}}')
    query_end_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($4==var2){print$6}else{print$3}}')
    echo -e $query"\t"$query_start_coords"\t"$query_end_coords"\t"$ref"\t"$ref_start_coords"\t"$ref_end_coords
    paf_syn_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window_${ref}/${ref}_${chr}.mapped-to.${query}_${chr}-${query_start_coords}_${query_end_coords}.syn.paf
    cat $paf_file | \
    awk -v var1=$query_start_coords -v var2=$query_end_coords '{if($8==var1&&$9==var2){print$0}}' > $paf_syn_file
    if [ -s $paf_syn_file ]
    then
    python /home/jl03308/git/NAM_pancentromere/NAM_SV/snp_dating4.py \
    -r $ref \
    -i1 $paf_syn_file \
    -i2 $aligned_syn_file \
    -s $query_start_coords -e $query_end_coords \
    -qs $ref_start_coords -qe $ref_end_coords \
    -o $output_file
    rm $paf_syn_file
    else
      echo "empty"
    fi
    done < $aligned_syn_file

fi

#done
fi
```

## Extract SNPs mapped to intergenic spaces

Find gene and UMR coordinates based on genome annotation
```
#!/bin/bash
chr="$1"
ref="$2"
ml BEDTools/2.29.2-GCC-8.3.0

NAMline=(B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
ref_lc=$(echo "$ref" | awk '{print tolower($0)}')
wdir=/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency
#chr=chr8
#ref=B73
chrlen=$(cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/$ref.*pseudomolecules-v*.fasta.fai | awk -v var1=$chr '{if($1==var1){print $2}}')
cd $wdir
cat ${chr}.allele_frequency.bed | \
awk -v var1=$ref -v OFS='\t' '{print var1,$0,$2-$1+1}' | \
sort -k1,1 -k4,4n | \
bedtools groupby -i - -g 1,4 -o sum -c 5 > ${chr}.allele_frequency.totallen.bed

#genes
cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/${ref_lc}/zea_mays${ref_lc}_core_3_87_1.gff | \
awk -v var1=$chr -v var2=$ref -v OFS='\t' '{if($3=="gene"&&$1==var1){print var2,$4,$5}}' | \
bedtools sort -i - | \
bedtools merge -i - > ${ref}.${chr}.gene.bed
#UMRs
cat /scratch/jl03308/NAM_pancentromere/methylation/UMRs/meth_${ref}.ref_${ref}.UMR.bed | \
awk '{print "chr"$0}' | \
awk -v var1=$chr -v var2=$ref -v OFS='\t' '{if($1==var1){print var2,$2,$3}}' | \
bedtools sort -i - | \
bedtools merge -i - > ${ref}.${chr}.UMR.bed

cat ${ref}.${chr}.gene.bed ${ref}.${chr}.UMR.bed |bedtools sort -i - | bedtools merge -i - > ${ref}.${chr}.UMR_gene.bed

totalgenelen=$(cat ${ref}.${chr}.gene.bed | awk '{sum+=$3-$2+1;} END{print sum;}')
totalUMRlen=$(cat ${ref}.${chr}.UMR.bed | awk '{sum+=$3-$2+1;} END{print sum;}')
totalgeneUMRlen=$(cat ${ref}.${chr}.UMR_gene.bed | awk '{sum+=$3-$2+1;} END{print sum;}')

```

Combine syntenic SNPs in vcf file called by `snp_dating.py` and remove SNPs mapped to gene and UMR regions
```
ml BEDTools/2.29.2-GCC-8.3.0
ref=B73
NAMline=(B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)

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

#map all syntenic snps to intergenic locations on each aligned blocks

for chr in chr{1..10}
do
  for line in ${NAMline[@]}
  do
    echo $line
    syn_snp_file=/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/$chr/${line}.${chr}.syn.snp.bed
    gene_UMR_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency/${ref}.${chr}.UMR_gene.bed
    syn_alignment_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/$chr/${ref}_${line}.syn.aligned.bed
  # remove syn snp from gene and UMR
  bedtools intersect -a $syn_snp_file -b $gene_UMR_file -v | \
  # map syn snps in intergenic space to syn aligned blocks
  bedtools intersect -a $syn_alignment_file -b - -c| \
  # calculate intergenic space in each aligned blocks
  bedtools intersect -a - -b $gene_UMR_file -wao | \
  bedtools sort -i - | \
  # print intergenic syn snps, and intergenic length in each aligned blocks
  # ref, start, end, query, start,end, intergenic snp number, gene+UMR length
  bedtools groupby -i - -g 1,2,3,4,5,6,14 -c 18 -o sum | \
  # ref, start, end, query, start,end, intergenic snp number, intergenic length
  awk -v OFS='\t' '{print $1,$2,$3,$4,$5,$6,$7,$3-$2-$8}' | \
  awk -v OFS='\t' -v var1=$chr '{if($8>0){print $0,$7/$8/2/3.3 *100000000,var1}}' >> /scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_div.bed
done
done
```

## Plot

```
ref='$1'

ml R
for chr in chr{1..10}
do
wdir=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}
ref_gff=/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${ref}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff
ref_genomesize=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${ref}.*pseudomolecules-v*.chrom.sizes
centromere_file=/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum
output_file=$wdir/${ref}_NAM.div.txt

if [ -f "$output_file" ]; then
    echo "$output_file exists."
    rm $output_file
fi

for file in $wdir/*${ref}*.aligned.div.txt
do
  cat $file | awk -v var1=$ref -v OFS='\t' '{if($1==var1){print$0}else{print $4,$5,$6,$1,$2,$3,$7,$8,$9}}'>> $output_file
done

Rscript /home/jl03308/git/NAM_pancentromere/NAM_SV/snp_dating_plot.R \
 $wdir \
 $ref_gff  \
 $centromere_file \
 $ref_genomesize \
 $output_file \
 ${chr} ${ref}
done
```
