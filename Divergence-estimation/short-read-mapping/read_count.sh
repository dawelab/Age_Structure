#!/bin/bash
#SBATCH --partition=highmem_p
#SBATCH --job-name=read_count
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=55:00:00
#SBATCH --mem=30gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jl03308@uga.edu

for file in /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/PC_*/*_2.fastq.gz
do
  prefix=$(basename $file | cut -f1 -d "_")
  dir=$(dirname $file)
  r1=$(zcat ${dir}/${prefix}_1.fastq.gz | echo $((`wc -l`/4)))
  r2=$(zcat $file | echo $((`wc -l`/4)))
  echo -e $prefix"\t"$r1"\t"$r2 >> /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/PC_readcount.txt
done

# count numbe of validated reads
for file in /scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/trimming/parviglumis/pc/SRR*.sh.e*
do
  prefix=$(basename $file |cut -f1 -d ".")
  val=$(cat $file |grep "Total number of sequences analysed:" |cut -f2 -d ":"| sed "s/^[ \t]*//")
  echo -e $prefix"\t"$val >> /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/PC_val_readcount.txt
done

for file in /scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/trimming/parviglumis/pc/pctrim_PC_*.sh.e*
do
  prefix=$(cat $file |grep "Input filename" | head -1 | cut -f2 -d ":"| cut -f1 -d "_" | sed "s/^[ \t]*//")
  val=$(cat $file |grep "Total number of sequences analysed:" |cut -f2 -d ":"| sed "s/^[ \t]*//")
  echo -e $prefix"\t"$val >> /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/PC_val_readcount.txt
done

cat /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/PC_val_readcount.txt |sort -k2 -n
