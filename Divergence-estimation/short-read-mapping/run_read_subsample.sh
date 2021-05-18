src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/subsample
mkdir -p $src_path
num=330000000
for PE1 in /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/PC_*/*_1_val_1.fq.gz
do
  prefix=$(basename $PE1 |cut -f1 -d "_")
  wdir=$(dirname $PE1)
  PE2=${wdir}/${prefix}_2_val_2.fq.gz
  file_name=$prefix.subsample.sh;
  echo '#!/bin/bash' >> ${src_path}/${file_name};
  echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
  echo '#SBATCH --job-name=subsample_'${prefix} >> ${src_path}/${file_name};
  echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
  echo '#SBATCH --ntasks=1'>> ${src_path}/${file_name};
  echo '#SBATCH --time=35:00:00'>> ${src_path}/${file_name};
  echo '#SBATCH --mem=100gb'>> ${src_path}/${file_name};
  echo '#SBATCH --mail-type=END,FAIL' >> ${src_path}/${file_name};
  echo '#SBATCH --mail-user=jl03308@uga.edu' >> ${src_path}/${file_name};
  echo 'sh ~/git/NAM_pancentromere/teosinte/read_subsample.sh' $PE1 $PE2 $num >> ${src_path}/${file_name};
done
