#!/bin/bash
chr="$1"
ref="$2"
query="$3"
#chr=chr1
#ref=B73
#query=P39
if [ $ref = $query ]; then
  echo "pass, 100% match"
else
module load BEDTools/2.29.2-GCC-8.3.0
module load Anaconda3/5.0.1
mkdir -p /scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window_${ref}
#allpaf=$(ls /scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/*_${chr}.mapped-to.${ref}_${chr}.sorted.paf |grep -v "Ab10")

# for paf_file in $allpaf
# do
#query=$(basename $paf_file |cut -f1 -d "." |cut -f1 -d "_")
FILE=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${ref}_${query}.syn.aligned.bed

if [ -f "$FILE" ]; then
    aligned_syn_file=$FILE
    paf_unsorted_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${query}_${chr}.mapped-to.${ref}_${chr}.paf
    paf_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${query}_${chr}.mapped-to.${ref}_${chr}.sorted.paf
    cat $paf_unsorted_file | sort -k6,6 -k8,8n > $paf_file
    #only include snps of syntenic aligned segments and inversions bigger than 20Kb
    output_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${ref}_${query}.aligned.div.txt
    if [ -f "$output_file" ]; then
        echo "$output_file exists."
        rm $output_file
    fi
    while IFS= read line
    do
    ref_start_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($1==var1){print$2}else{print$5}}')
    ref_end_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($1==var1){print$3}else{print$6}}')
    query_start_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($4==var2){print$5}else{print$2}}')
    query_end_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($4==var2){print$6}else{print$3}}')
    echo -e $ref"\t"$ref_start_coords"\t"$ref_end_coords"\t"$query"\t"$query_start_coords"\t"$query_end_coords
    paf_syn_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window_${ref}/${query}_${chr}.mapped-to.${ref}_${chr}-${ref_start_coords}_${ref_end_coords}.syn.paf
    cat $paf_file | \
    awk -v var1=$ref_start_coords -v var2=$ref_end_coords '{if($8==var1&&$9==var2){print$0}}' > $paf_syn_file
    if [ -s $paf_syn_file ]
    then
    python /home/jl03308/git/NAM_pancentromere/NAM_SV/snp_dating4.py \
    -r $ref \
    -i1 $paf_syn_file \
    -i2 $aligned_syn_file \
    -s $ref_start_coords -e $ref_end_coords \
    -qs $query_start_coords -qe $query_end_coords \
    -o $output_file
    rm $paf_syn_file
    else
      echo "empty"
    fi
    done < $aligned_syn_file

else
    aligned_syn_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${query}_${ref}.syn.aligned.bed
    paf_unsorted_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${ref}_${chr}.mapped-to.${query}_${chr}.paf
    paf_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/${ref}_${chr}.mapped-to.${query}_${chr}.sorted.paf
    cat $paf_unsorted_file | sort -k6,6 -k8,8n > $paf_file
    #only include snps of syntenic aligned segments and inversions bigger than 20Kb
    output_file=/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/${chr}/${query}_${ref}.aligned.div.txt
    if [ -f "$output_file" ]; then
        echo "$output_file exists."
        rm $output_file
    fi
    while IFS= read line
    do
    ref_start_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($1==var1){print$2}else{print$5}}')
    ref_end_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($1==var1){print$3}else{print$6}}')
    query_start_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($4==var2){print$5}else{print$2}}')
    query_end_coords=$(echo $line |awk -v var1=$ref -v var2=$query '{if($4==var2){print$6}else{print$3}}')
    echo -e $query"\t"$query_start_coords"\t"$query_end_coords"\t"$ref"\t"$ref_start_coords"\t"$ref_end_coords
    paf_syn_file=/scratch/jl03308/NAM_pancentromere/genome_alignment/${chr}_shujun/snp_window_${ref}/${ref}_${chr}.mapped-to.${query}_${chr}-${query_start_coords}_${query_end_coords}.syn.paf
    cat $paf_file | \
    awk -v var1=$query_start_coords -v var2=$query_end_coords '{if($8==var1&&$9==var2){print$0}}' > $paf_syn_file
    if [ -s $paf_syn_file ]
    then
    python /home/jl03308/git/NAM_pancentromere/NAM_SV/snp_dating4.py \
    -r $ref \
    -i1 $paf_syn_file \
    -i2 $aligned_syn_file \
    -s $query_start_coords -e $query_end_coords \
    -qs $ref_start_coords -qe $ref_end_coords \
    -o $output_file
    rm $paf_syn_file
    else
      echo "empty"
    fi
    done < $aligned_syn_file

fi

#done
fi
