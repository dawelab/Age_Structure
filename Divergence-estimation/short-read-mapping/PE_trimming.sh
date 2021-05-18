#!/bin/bash
PE1="$1"
PE2="$2"

wdir=$(dirname $PE1)

module load Trim_Galore/0.6.5-GCCcore-8.3.0-Java-11-Python-3.7.4
cd $wdir
trim_galore --fastqc --gzip --paired $PE1 $PE2 -o $wdir
