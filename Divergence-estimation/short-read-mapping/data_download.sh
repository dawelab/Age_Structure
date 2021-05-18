#!/bin/bash
sra="$1"
out_dir="$2"

module load SRA-Toolkit/2.9.6-1-centos_linux64
fasterq-dump $sra -e 12 -f -S -p -O $out_dir
fastq_files=$(find $out_dir -name '*.fastq')
for file in $fastq_files
do
  gzip $file
done
