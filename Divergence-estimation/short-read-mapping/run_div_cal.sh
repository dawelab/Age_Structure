src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/div_cal
mkdir -p $src_path
genomecovfiles=$(ls /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/*/*.effective.genomecov)
for file in $genomecovfiles
do
  prefix=$(basename $file |cut -f1 -d ".")
  wdir=$(dirname $file)
  #bedtools merge -i $file > ${wdir}/${prefix}.mask.bed
  file_name=$prefix.div_cal.sh;
  snp_file=$(ls ${wdir}/*.flt.snps.bed |grep -v "P39")
  echo '#!/bin/bash' >> ${src_path}/${file_name};
  echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
  echo '#SBATCH --job-name='${prefix}'_bedtoolscov' >> ${src_path}/${file_name};
  echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
  echo '#SBATCH --ntasks=2'>> ${src_path}/${file_name};
  echo '#SBATCH --time=10:00:00'>> ${src_path}/${file_name};
  echo '#SBATCH --mem=20gb'>> ${src_path}/${file_name};
  echo '#SBATCH --mail-type=END,FAIL' >> ${src_path}/${file_name};
  echo '#SBATCH --mail-user=jl03308@uga.edu' >> ${src_path}/${file_name};
  echo 'sh ~/git/NAM_pancentromere/teosinte/div_cal.sh' $file $snp_file $ref >> ${src_path}/${file_name};
done

#dir=$(find /scratch/jl03308/NAM_pancentromere/rawdata/teosinte/parviglumis/ -type d -iname "TIL*")
#bamfiles=$(ls /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/*/*.q20.sorted.bam |grep -v "PC" | grep -v "P39")
src_path=/scratch/jl03308/NAM_pancentromere/src/run_scripts/teosinte/div_cal
ref=P39
genomecovfiles=$(ls /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/*/*.P39.effective.genomecov)
for file in $genomecovfiles
do
  prefix=$(basename $file |cut -f1 -d ".")
  wdir=$(dirname $file)
  #bedtools merge -i $file > ${wdir}/${prefix}.mask.bed
  file_name=$prefix.P39.div_cal.sh;
  snp_file=$(ls ${wdir}/*.P39.flt.snps.bed)
  echo '#!/bin/bash' >> ${src_path}/${file_name};
  echo '#SBATCH --partition=batch' >> ${src_path}/${file_name};
  echo '#SBATCH --job-name='${prefix}'_bedtoolscov' >> ${src_path}/${file_name};
  echo '#SBATCH --nodes=1'>> ${src_path}/${file_name};
  echo '#SBATCH --ntasks=2'>> ${src_path}/${file_name};
  echo '#SBATCH --time=10:00:00'>> ${src_path}/${file_name};
  echo '#SBATCH --mem=20gb'>> ${src_path}/${file_name};
  echo '#SBATCH --mail-type=END,FAIL' >> ${src_path}/${file_name};
  echo '#SBATCH --mail-user=jl03308@uga.edu' >> ${src_path}/${file_name};
  echo 'sh ~/git/NAM_pancentromere/teosinte/div_cal.sh' $file $snp_file $ref >> ${src_path}/${file_name};
done
