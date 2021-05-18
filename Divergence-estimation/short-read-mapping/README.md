## Tools used for SNP calling and divergence dating
- Read processing
    * subsampling - seqtk (v1.3)
    * trimming and QC - TrimGalore (v0.6.5)
    * mapping - BWA (v0.7.17)
    * filtering - samtools (v1.9)
- SNP calling
    * bcftools mpileup (v1.6)
- Divergence calculation and normalization
    * snpable region - bedtools (2.29.2)
- Plot
    * karyoploteR

## Read processing (paired-end)
subsampling: `read_subsample.sh`
trimming: `PE_trimming.sh`
mapping and filtering: `PE_alignment.sh`

## SNP calling
`mpileup.sh` performs SNP and small indel calling with BCFtools mpileup, and filtering based on confidence.

## Divergence calculation and normalization

Step1: `bedtoolscov.sh` performs coverage and SNPable area calculation in a specific window, and outputs `$prefix.effective.genomecov`, which is used for divergence calculation in the next step.
Example:
```
bamfiles=$(ls /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/*/*.P39.q20.sorted.bam)
for file in bamfiles
do
sh bedtoolscov.sh' $file
done
```

Step2: `div_cal.sh` runs coverage and SNPable area calculation in a specific window,
snp extraction in intergenic areas and divergence time calculation.
Example:
```
ref=P39
genomecovfiles=$(ls /scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/*/*.P39.effective.genomecov)
for file in $genomecovfiles
do
  prefix=$(basename $file |cut -f1 -d ".")
  wdir=$(dirname $file)
  #bedtools merge -i $file > ${wdir}/${prefix}.mask.bed
  file_name=$prefix.P39.div_cal.sh;
  snp_file=$(ls ${wdir}/*.P39.flt.snps.bed)
  sh div_cal.sh' $file $snp_file $ref
done
```

## Plot

```
Rscript snp_plot.normalized.R \
$wdir \
$ref_gff  \
$centromere_file \
$ref_genomesize \
$output_file \
${ref} \
${chr} \
Teo_NAM_${ref}
```
