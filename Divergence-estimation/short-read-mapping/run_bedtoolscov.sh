src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/bedtools_cov
mkdir -p $src_path
bamfiles=$(ls /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/PC*/*.q20.sorted.bam |grep -v "P39")
for file in $bamfiles
do
  prefix=$(basename $file |cut -f1 -d ".")
  wdir=$(dirname $file)
  file_name=$prefix.bedtools_cov.sh;
  echo '#!/bin/bash' >> ${src_path}/${file_name};
  echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
  echo '#SBATCH --job-name='${prefix}'_bedtoolscov' >> ${src_path}/${file_name};
  echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
  echo '#SBATCH --ntasks=2'>> ${src_path}/${file_name};
  echo '#SBATCH --time=10:00:00'>> ${src_path}/${file_name};
  echo '#SBATCH --mem=20gb'>> ${src_path}/${file_name};
  echo '#SBATCH --mail-type=END,FAIL' >> ${src_path}/${file_name};
  echo '#SBATCH --mail-user=jl03308@uga.edu' >> ${src_path}/${file_name};
  echo 'sh ~/git/NAM_pancentromere/teosinte/bedtoolscov.sh' $file >> ${src_path}/${file_name};
done

#dir=$(find /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/ -type d -iname "TIL*")

src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/bedtools_cov
mkdir -p $src_path
bamfiles=$(ls /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/*/*.P39.q20.sorted.bam | grep -v "TIL" |grep -v "PC")
for file in $bamfiles
do
  prefix=$(basename $file |cut -f1 -d ".")
  wdir=$(dirname $file)
  file_name=$prefix.P39.bedtools_cov.sh;
  echo '#!/bin/bash' >> ${src_path}/${file_name};
  echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
  echo '#SBATCH --job-name='${prefix}'_bedtoolscov' >> ${src_path}/${file_name};
  echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
  echo '#SBATCH --ntasks=2'>> ${src_path}/${file_name};
  echo '#SBATCH --time=10:00:00'>> ${src_path}/${file_name};
  echo '#SBATCH --mem=20gb'>> ${src_path}/${file_name};
  echo '#SBATCH --mail-type=END,FAIL' >> ${src_path}/${file_name};
  echo '#SBATCH --mail-user=jl03308@uga.edu' >> ${src_path}/${file_name};
  echo 'sh ~/git/NAM_pancentromere/teosinte/bedtoolscov.sh' $file >> ${src_path}/${file_name};
done
