src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/read_alignment/parviglumis/
mkdir -p $src_path
for R1 in /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/PC_*/*_1_val_1.330000000.fastq.gz
 do
   prefix=$(basename $R1 | cut -f1 -d "_")
   dir=$(dirname $R1)
   echo $prefix
   R2=${dir}/${prefix}_2_val_2.330000000.fastq.gz
   # if [ -f "${trimmed}" ]
   # then
   #   echo 'Raw reads already trimmed'
  #else
    file_name=$prefix.sh;
    echo '#!/bin/bash' >> ${src_path}/${file_name};
    echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
    echo '#SBATCH --job-name='${prefix}'_PEalignment' >> ${src_path}/${file_name};
    echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
    echo '#SBATCH --ntasks=12'>> ${src_path}/${file_name};
    echo '#SBATCH --time=35:00:00'>> ${src_path}/${file_name};
    echo '#SBATCH --mem=100gb'>> ${src_path}/${file_name};
    echo '#SBATCH --mail-type=END,FAIL' >> ${src_path}/${file_name};
    echo '#SBATCH --mail-user=jl03308@uga.edu' >> ${src_path}/${file_name};
    echo 'sh /home/jl03308/git/NAM_pancentromere/teosinte/PE_alignment.sh' $R1 $R2 >> ${src_path}/${file_name};
#  fi
done
