#!/bin/bash

R1="$1"
R2="$2"
ref_genome="$3"

prefix=$(basename $R1 |cut -f1 -d "_")
ref=$(basename $ref_genome |cut -f1 -d ".")

module load BWA/0.7.17-GCC-8.3.0
mkdir -p /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${prefix}
cd /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${prefix}

index=$ref_genome.ann
if [ -f "${index}" ]; then
   echo "skip genome indexing."
else
  bwa index $ref_genome
fi

bwa mem -t 12 $ref_genome $R1 $R2 > ${prefix}.${ref}.sam

module load BWA/0.7.17-GCC-8.3.0
mkdir -p /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${prefix}
cd /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/${prefix}

module load SAMtools/1.9-GCC-8.3.0
samtools view -@ 12 -b -o ${prefix}.${ref}.bam ${prefix}.${ref}.sam
samtools sort -o ${prefix}.${ref}.sorted.bam -T ${prefix}.${ref} -@ 12 ${prefix}.${ref}.bam
samtools index -@ 12 ${prefix}.${ref}.sorted.bam
samtools view -@ 12 -bhq 20 ${prefix}.${ref}.sorted.bam -o ${prefix}.${ref}.q20.bam
samtools sort -o ${prefix}.${ref}.q20.sorted.bam -T ${prefix}.${ref} -@ 12 ${prefix}.${ref}.q20.bam
samtools index -@ 12 ${prefix}.${ref}.q20.sorted.bam

# module load picard/2.21.6-Java-11

# module load picard/2.21.6-Java-11
# java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=${prefix}.sorted.bam O=${prefix}.rmdup.bam M=${prefix}.marked_dup_metrics.txt REMOVE_DUPLICATES=true
# samtools sort -o ${prefix}.rmdup.sorted.bam -T ${prefix} -@ 12 ${prefix}.rmdup.bam
# samtools index -@ 12 ${prefix}.rmdup.sorted.bam
# samtools view -@ 12 -bhq 20 ${prefix}.rmdup.sorted.bam -o ${prefix}.rmdup.q20.bam
# samtools sort -o ${prefix}.rmdup.q20.sorted.bam -T ${prefix} -@ 12 ${prefix}.rmdup.q20.bam
# samtools flagstat ${prefix}.rmdup.q20.sorted.bam > ${prefix}.flagstat
# samtools index -@ 12 ${prefix}.rmdup.q20.sorted.bam
rm ${prefix}.${ref}.sam
rm ${prefix}.${ref}.bam
# rm ${prefix}.rmdup.bam
# rm ${prefix}.rmdup.q20.bam
