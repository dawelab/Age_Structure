#!/bin/bash
line="$1"
chr="$2"
start="$3"
end="$4"
type="$5"
prefix="$6"

#start="$3"
#end="$4"
if [[ $type == "centc" ]]
then
  repeat="Cent/CentC"
  #echo $repeat
elif [[ $type == "knob180" ]]; then
  repeat="knob/knob180"
  #echo $repeat
elif [[ $type == "TR-1" ]]; then
  repeat="knob/TR-1"
  #echo $repeat
elif [[ $type == "rDNA" ]]; then
  repeat="rDNA/spacer"
  #echo $repeat
elif [[ $type == "subtelomere" ]]; then
  repeat="subtelomere/4-12-1"
  #echo $repeat
fi

module load SAMtools/1.9-GCC-8.3.0
module load BEDTools/2.29.2-GCC-8.3.0

mkdir -p /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}
cd /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}

repeat_bed=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${line}.bed
repeat_fa=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${line}.fa

#extract index each repeat in selected region
cat /scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/${line}.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff | \
awk -v var1=$repeat '{if($3==var1){print$0}}' | \
awk -v var3=${line}_${chr} -v var1=$start -v var2=$end '{if($1==var3&&$4>=var1&&$5<=var2){print$0}}'| \
awk -v var2=$chr '{print var2"\t"$4"\t"$5"\t"$7}' |cat -n -| \
awk -v var1=$type '{print $2"\t"$3"\t"$4"\t"var1"."$1"\t0\t"$5}' > $repeat_bed
#extract repeat sequence
bedtools getfasta -fi /scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/${line}.*pseudomolecules-v*.fasta -bed $repeat_bed -nameOnly -s > $repeat_fa
sed -i 's/(.*)//g' $repeat_fa
if [ -s "$repeat_fa" ]
then
samtools faidx $repeat_fa
fi
#split repeat seq file into individual monomer files
mkdir -p /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${prefix}/${line}
cd /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${prefix}/${line}
if [ -s "$repeat_fa" ]
then
repeat_index=$(cat $repeat_fa |grep ">" |sed 's|>||g')
for x in $repeat_index
do
samtools faidx $repeat_fa $x > ${type}.${chr}.${prefix}.${line}.$x.fa
done
fi
