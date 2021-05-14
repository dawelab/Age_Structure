wdir=/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating
cd $wdir
knobs=$(ls NAM.*.chr*.knob*.pairwise_distance.txt | cut -f3,4 -d "." |sort |uniq)
for prefix in $knobs
do
file1=NAM.knob180.${prefix}.pairwise_distance.txt
file2=NAM.TR-1.${prefix}.pairwise_distance.txt
out=NAM.knob180+TR-1.${prefix}.pairwise_distance.txt
if [ -f "$file1" ] && [ -f "$file2" ]; then
#join -j 4 -o 1.4,1.5,1.6,1.7,1.8,1.9,2.6,2.7,2.8,2.9 $file1 $file2 | head
Rscript /home/jl03308/git/NAM_pancentromere/syn_repeat/merge_files.R $file1 $file2 $out
Rscript ~/git/NAM_pancentromere/syn_repeat/knob_cluster.R  \
$wdir \
$out \
NAM.knob180+TR-1.${prefix} 
fi
done
