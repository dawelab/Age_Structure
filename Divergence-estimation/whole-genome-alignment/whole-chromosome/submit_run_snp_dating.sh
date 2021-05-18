ref=CML247
ref=Ki3

src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/NAM_SVs/NAM_${ref}
mkdir $src_path
cd $src_path
lines=($(cat /scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centC.array.assembly.status.content | cut -f1 |sort |uniq|grep -v "B73" |grep -v $ref))
for query in "${lines[@]}"
do
  for chr in chr{1..10}
  do

  mkdir -p $src_path
  file_name=${query}.${chr}.${ref}.sh;
  echo '#!/bin/bash' >> ${src_path}/${file_name};
  echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
  echo '#SBATCH --job-name='${query}.${chr}.${ref} >> ${src_path}/${file_name};
  echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
  echo '#SBATCH --ntasks=1'>> ${src_path}/${file_name};
  echo '#SBATCH --time=35:00:00'>> ${src_path}/${file_name};
  echo '#SBATCH --mem=30gb'>> ${src_path}/${file_name};
  echo 'sh ~/git/NAM_pancentromere/NAM_SV/run_snp_dating5.sh' $chr $ref $query >> ${src_path}/${file_name};
done
done
