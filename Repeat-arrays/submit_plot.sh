type=TR-1
type=knob180

cd /lustre2/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/arrays/
chrs=$(ls *knob*.bed | cut -f1 -d "."|sort |uniq)
for chr in $chrs
do
  echo $chr
  src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/syn_repeat/plot/${chr}
  mkdir -p $src_path
  prefixes=$(ls ${chr}.knob*.bed | cut -f2 -d ".")
  for prefix in $prefixes
  do
  echo $prefix
  file_name=${type}.${chr}.${prefix}.sh;
  echo '#!/bin/bash' >> ${src_path}/${file_name};
  echo '#SBATCH --partition=highmem_p' >> ${src_path}/${file_name};
  echo '#SBATCH --job-name='${type}.${chr}.${prefix} >> ${src_path}/${file_name};
  echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
  echo '#SBATCH --ntasks=1'>> ${src_path}/${file_name};
  echo '#SBATCH --time=55:00:00'>> ${src_path}/${file_name};
  echo '#SBATCH --mem=140gb'>> ${src_path}/${file_name};
  echo 'sh ~/git/NAM_pancentromere/syn_repeat/plot_pairwise_repeats.sh' $type $chr $prefix >> ${src_path}/${file_name};
done
done

type=centc

cd /lustre2/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/arrays/
chrs=$(ls *centc*.bed | cut -f1 -d "."|sort |uniq)
for chr in $chrs
do
  echo $chr
  src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/syn_repeat/plot/${chr}
  mkdir -p $src_path
  prefixes=$(ls ${chr}.centc*.bed | cut -f2 -d ".")
  for prefix in $prefixes
  do
  echo $prefix
  file_name=${type}.${chr}.${prefix}.sh;
  echo '#!/bin/bash' >> ${src_path}/${file_name};
  echo '#SBATCH --partition=highmem_p' >> ${src_path}/${file_name};
  echo '#SBATCH --job-name='${type}.${chr}.${prefix} >> ${src_path}/${file_name};
  echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
  echo '#SBATCH --ntasks=1'>> ${src_path}/${file_name};
  echo '#SBATCH --time=55:00:00'>> ${src_path}/${file_name};
  echo '#SBATCH --mem=140gb'>> ${src_path}/${file_name};
  echo 'sh ~/git/NAM_pancentromere/syn_repeat/plot_pairwise_repeats.sh' $type $chr $prefix >> ${src_path}/${file_name};
done
done
