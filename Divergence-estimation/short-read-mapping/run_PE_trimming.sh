src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/trimming/parviglumis/pc
mkdir -p $src_path
for x in /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/PC_*/*_2.fastq.gz
 do
   prefix=$(basename $x | cut -f1 -d "_")
   dir=$(dirname $x)
   echo $prefix
   trimmed=${dir}/${prefix}_1_val_1.fq.gz
   if [ -f "${trimmed}" ]
   then
     echo 'Raw reads already trimmed'
  else
    file_name=$prefix.sh;
    echo '#!/bin/bash' >> ${src_path}/${file_name};
    echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
    echo '#SBATCH --job-name=PEtrim_'${prefix} >> ${src_path}/${file_name};
    echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
    echo '#SBATCH --ntasks=2'>> ${src_path}/${file_name};
    echo '#SBATCH --time=35:00:00'>> ${src_path}/${file_name};
    echo '#SBATCH --mem=60gb'>> ${src_path}/${file_name};
    echo '#SBATCH --mail-type=END,FAIL' >> ${src_path}/${file_name};
    echo '#SBATCH --mail-user=jl03308@uga.edu' >> ${src_path}/${file_name};
    echo 'sh /scratch/jl03308/NAM_pancentromere/src/core_scripts/PE_trimming.sh' ${dir}/${prefix}_1.fastq.gz ${dir}/${prefix}_2.fastq.gz >> ${src_path}/${file_name};
  fi
done

dir=$(find /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/ -type d -iname "TIL*")

for x in $dir
do
dirname=$(basename $x)
cat $x/*_1.fastq.gz > $x/${dirname}_1.fastq.gz
cat $x/*_2.fastq.gz > $x/${dirname}_2.fastq.gz
mkdir -p $x/SRA
mv $x/SRR* $x/SRA
done
src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/trimming/parviglumis/
for x in /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/Mex_ISU*/*_1.fastq.gz
 do
   prefix=$(basename $x | cut -f1 -d "_")
   dir=$(dirname $x)
   echo $prefix
   trimmed=${dir}/${prefix}_1_val_1.fq.gz
   if [ -f "${trimmed}" ]
   then
     echo 'Raw reads already trimmed'
  else
    file_name=$prefix.sh;
    echo '#!/bin/bash' >> ${src_path}/${file_name};
    echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
    echo '#SBATCH --job-name=PEtrim_'${prefix} >> ${src_path}/${file_name};
    echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
    echo '#SBATCH --ntasks=2'>> ${src_path}/${file_name};
    echo '#SBATCH --time=35:00:00'>> ${src_path}/${file_name};
    echo '#SBATCH --mem=60gb'>> ${src_path}/${file_name};
    echo '#SBATCH --mail-type=END,FAIL' >> ${src_path}/${file_name};
    echo '#SBATCH --mail-user=jl03308@uga.edu' >> ${src_path}/${file_name};
    echo 'sh /scratch/jl03308/NAM_pancentromere/src/core_scripts/PE_trimming.sh' ${dir}/${prefix}_1.fastq.gz ${dir}/${prefix}_2.fastq.gz >> ${src_path}/${file_name};
  fi
done
