# identify alignment range for figure 1 (46Mb to 51Mb) 
start=46000000 
end=51000000 
cat NAM.chr8.pairwise.syn.bed | \
awk -v var1=$start -v var2=$end '{if($2>var1&&$3<=var2&&$1=="B73"&&$4=="CML322"){print$0}}' > B73_CML322.syn.selected.bed

# redefine start and end coords
CML322_start=$(cat B73_CML322.syn.selected.bed| sort -k2,2n |cut -f2 | head -1)
CML322_end=$(cat B73_CML322.syn.selected.bed| sort -k2,2n |cut -f3 | tail -1)
echo -e "B73\t"$CML322_start"\t"$CML322_end > B73_CML322.range.bed
# find genes specific to B73
bedtools subtract -a B73_CML322.range.bed -b B73_CML322.syn.selected.bed > B73_CML322.B73specific.bed

# reformat the TE insertion time file, and find insertion time belonging to file 'B73_CML322.B73specific.bed'
# output format : B73chr8, start coord, end coord, SuperFamily, TE type, insertion time
cat B73.PLATINUM.pseudomolecules-v1.fasta.mod.pass.list | sed 's|B73_||g' | \
awk '{if($1~"chr8"){print$0}}' | cut -f1,10,11,12 |sed 's|chr8:|B73\t|g' | \
sed 's/[.]/\t/g' | \
cut -f1,2,4,5,6,7,8 | \
bedtools intersect -a - -b B73_CML322.B73specific.bed -wa > B73_CML322.B73specific.insertiontime.bed
