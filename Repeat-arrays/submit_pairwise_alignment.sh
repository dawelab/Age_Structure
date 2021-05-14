
for syn_repeat_file in /scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/arrays/chr*.bed
 do
   chr=$(basename $syn_repeat_file | cut -f1 -d ".")
   prefix=$(basename $syn_repeat_file | cut -f2 -d ".")
   repeat_type=$(basename $syn_repeat_file | cut -f2 -d "." |sed 's/[0-9]*//g')
   src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/syn_repeat/pairwise/${chr}/
   mkdir -p $src_path
   if [[ $repeat_type == "knob" ]]
   then
     types=(knob180 TR-1)
   else
     types=($repeat_type)
   fi

   for type in "${types[@]}"
   do
     echo $type
   #pairwise-alignment of knobs
   all_line=($(cat $syn_repeat_file |cut -f1 |sort |uniq))
   len=${#all_line[@]}
   ## Use bash for loop
   for ((i = 0; i < len; i++)); do
     for ((j = i+1; j < len; j++)); do
       line1=${all_line[$i]}
       line2=${all_line[$j]}
       #echo "$line1,$line2"
       input_fa1=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${line1}.fa
       input_fa2=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/seq/${type}/${chr}/${type}.${chr}.${prefix}.${line2}.fa
       if [ -s $input_fa1 ] && [ -s $input_fa2 ]
       then
       file_name=${type}.${chr}.${prefix}.${line1}_${line2}.sh;
       echo '#!/bin/bash' >> ${src_path}/${file_name};
       echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
       echo '#SBATCH --job-name='${line1}_${line2}.${type}.${chr}.${prefix} >> ${src_path}/${file_name};
       echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
       if [[ $type == "rDNA" ]]
       then
       echo '#SBATCH --ntasks=20'>> ${src_path}/${file_name};
       echo '#SBATCH --time=105:00:00'>> ${src_path}/${file_name};
       echo '#SBATCH --mem=1gb'>> ${src_path}/${file_name};
       echo 'sh /home/jl03308/git/NAM_pancentromere/syn_repeat/pairwise_alignment.sh' $input_fa1 $input_fa2 >> ${src_path}/${file_name};
      else
        echo '#SBATCH --ntasks=1'>> ${src_path}/${file_name};
        echo '#SBATCH --time=40:00:00'>> ${src_path}/${file_name};
        echo '#SBATCH --mem=1gb'>> ${src_path}/${file_name};
        echo 'sh /home/jl03308/git/NAM_pancentromere/syn_repeat/pairwise_alignment.sh' $input_fa1 $input_fa2 "BLAT">> ${src_path}/${file_name};
     fi
   fi
     done
   done
   done

done
