src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/deeptools
mkdir -p $src_path
for file in /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/PC*/*.q20.sorted.bam
do
  prefix=$(basename $file |cut -f1 -d ".")
  wdir=$(dirname $file)
  file_name=$prefix.deeptools.sh;
  echo '#!/bin/bash' >> ${src_path}/${file_name};
  echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
  echo '#SBATCH --job-name=bamcal_'${prefix} >> ${src_path}/${file_name};
  echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
  echo '#SBATCH --ntasks=12'>> ${src_path}/${file_name};
  echo '#SBATCH --time=15:00:00'>> ${src_path}/${file_name};
  echo '#SBATCH --mem=30gb'>> ${src_path}/${file_name};
  echo '#SBATCH --mail-type=END,FAIL' >> ${src_path}/${file_name};
  echo '#SBATCH --mail-user=jl03308@uga.edu' >> ${src_path}/${file_name};
  echo 'sh ~/git/NAM_pancentromere/teosinte/deeptools.sh' $file >> ${src_path}/${file_name};
done
