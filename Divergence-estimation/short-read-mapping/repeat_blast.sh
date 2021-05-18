
module load BLAST/2.2.26-Linux_x86_64

cd /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats

genome=/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/diplo.asm.p_chr_purged.fasta

formatdb -p F -i $genome

blastall -p blastn -d $genome -i /home/jl03308/reference/knob180.fa -m 8 > knob180.blast
blastall -p blastn -d $genome -i /home/jl03308/reference/CentC.fa -m 8 > CentC.blast
blastall -p blastn -d $genome -i /home/jl03308/reference/TR-1.fa -m 8 > TR-1.blast
blastall -p blastn -d $genome -i /home/jl03308/reference/rDNA_intergenic_spacer.fa -m 8 > rDNA.blast
blastall -p blastn -d $genome -i /home/jl03308/reference/subtelomeric_4-12-1.fa -m 8 > subtelomere.blast
blastall -p blastn -d $genome -i /home/jl03308/reference/cenH3.fa -m 8 > cenH3.blast

line1=$(basename $genome | cut -f1 -d ".")

chr=chr3
cat /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC.blast | grep $chr | sort -k9,9n | awk -v OFS='\t' '{if($4>30){print$2,$9,$10}}'| \
awk -v OFS='\t' '{if($3>$2){print$1,$2,$3,"+"}else{print$1,$3,$2,"-"}}' | awk '{print $0"\t"NR}' > \
/scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/CentC.${chr}.bed

module load BEDTools/2.29.2-GCC-8.3.0

'''extract centc, label centc by line, chr, and index, concatenate into one line, multiple sequence alignment, and phylogenetic treee'''
cat /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/CentC.${chr}.bed | \
awk -v OFS='\t' -v var2=$chr '{print $1,$2,$3,"diplo_"var2"_centc"$5,1,$4}' | \
bedtools getfasta -fi $ref -bed - -nameOnly -s > \
/scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/CentC.${chr}.named.fa

ml SAMtools/1.9-GCC-8.3.0
index=$(cat /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/CentC.${chr}.named.fa| grep ">" | sed 's|>||g')
mkdir -p /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/${chr}/seq
cd /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/${chr}/seq
for x in $index
do
  prefix=$(echo $x | cut -f1 -d "(")
samtools faidx /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/CentC.${chr}.named.fa $x > diplo.${chr}.${prefix}.fa
done

ml BLAT/3.5-GCC-8.3.0
mkdir -p /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/${chr}/diplo_${line}
line_repeat_fasta=/scratch/jl03308/NAM_pancentromere/analysis/centc_structure/identity_analysis_blat/${line}/${chr}/${line}.${chr}.centc.named.fa
for x in $index
do
  prefix=$(echo $x | cut -f1 -d "(")
  blat -minIdentity=70 -maxGap=10 -minScore=0 -repMatch=2147483647 diplo.${chr}.${prefix}.fa $line_repeat_fasta /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/${chr}/diplo_${line}/${line}.${chr}.${prefix}.psl
done

cat /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/${chr}/diplo_${line}/${line}.${chr}.*.psl | \
grep "centc" | cut -f1,10,11,14,15  | sed 's|_|\t|g' | cut -f1,4,5,8,9 | \
sed 's|centc||g' |sed 's|(+)||g' | sed 's|(-)||g' | \
awk '{if($2>$4){print$1"\t"$4"\t"$5"\t"$2"\t"$3}else{print$0}}' | \
awk '{if($2<$4){print$2"\t"$4"\t"$1*100/($3+$5-$1)}}' | \
sort -k1,1n -k2,2n -u > /scratch/jl03308/NAM_pancentromere/teosinte/diploperennis/repeats/CentC/${chr}/diplo_${line}.${chr}.centc.stat.plot
#cat ${line}.${chr}.centc.psl | grep "centc" | cut -f1,10,11,14,15  |sed 's|centc||g' |awk '{if($2<$4){print$2"\t"$4"\t"$1*100/($3+$5-$1)}}' |sort -k1,1n -k2,2n > ${line}.${chr}.centc.stat.plot
# pairwise alignment identity summary for distribution plot (violin + boxplot)
#cat ${line}.${chr}.centc.stat.plot | awk -v var1=$line -v var2=$chr '{print var1"\t"var2"\t"$3}' > ${line}.${chr}.stats
# pairwise alignment identity summary for distribution plot (violin + boxplot)
cat ${line}.${chr}.centc.stat.plot | awk -v var1=$line -v var2=$chr '{print var1"\t"var2"\t"$0}' > ${line}.${chr}.labelled.stat.plot
