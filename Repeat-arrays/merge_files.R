#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
file1 <- args[1]
file2 <- args[2]
out <- args[3]
#file1<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/NAM.knob180.chr5.knob1.pairwise_distance.txt'
#file2<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/NAM.TR-1.chr5.knob1.pairwise_distance.txt'
#out<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/NAM.knob180+TR-1.chr5.knob1.pairwise_distance.txt'

if (file.exists(file1))
  file1<-read.table(file1)
if (file.exists(file2))
  file2<-read.table(file2)

file_merged <- merge(file1,file2,by=c("V4","V5"))
file_merged$reflen <- file_merged$V6.x + file_merged$V6.y
file_merged$querylen <- file_merged$V7.x + file_merged$V7.y
file_merged$totalsimilarity <- file_merged$V8.x + file_merged$V8.y
file_merged$totalalignlen <- file_merged$V9.x + file_merged$V9.y
file_merged$meansimilarity <- file_merged$totalsimilarity/file_merged$totalalignlen
file_merged$V1 <- unique(file1$V1)
file_merged$V2 <- unique(file1$V2)
file_merged$V3 <- 'knob180+TR-1'

file_merged <- file_merged[c('V1','V2','V3','V4','V5','reflen','querylen','totalsimilarity','totalalignlen','meansimilarity')]
write.table(file_merged, file = out, sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
