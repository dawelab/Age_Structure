output=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/NAM.repeat.num.txt
NAMline=(B73 B97 CML103 CML228 CML247 CML277 CML322 CML333 CML52 CML69 HP301 IL14H Ki11 Ki3 Ky21 M162W M37W Mo18W MS71 NC350 NC358 Oh43 Oh7b P39 Tx303 Tzi8)
repeattype=(centc knob180 TR-1 rDNA)
for type in ${repeattype[@]}
do
  chrs=$(ls /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/*/*.fa |cut -f2 -d "." |sort |uniq)
for chr in $chrs
do
repeatprefix=$(ls /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.*.fa |cut -f3 -d "." |sort |uniq)
for prefix in $repeatprefix
do
for line in ${NAMline[@]}
  do
  seq_file=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${line}.fa
  if [ -f "$seq_file" ]; then
  num=$(cat $seq_file |grep -c ">")
  echo -e $type"\t"$chr"\t"$prefix"\t"$line"\t"$num >> $output
  else
  echo -e $type"\t"$chr"\t"$prefix"\t"$line"\t0" >> $output
  fi
  done
done
done
done
