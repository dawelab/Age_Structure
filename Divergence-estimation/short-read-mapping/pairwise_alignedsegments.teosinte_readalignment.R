#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
wd_dir <- args[1]
paf_file <- args[2]
centromere_file <- args[3]
gff_file1 <- args[4]
gff_file2 <- args[5]
gap_file <- args[6]
ref<-args[7]
query<-args[8]
chr <- args[9]
alignment_file <- args[10]

# wd_dir <- '/scratch/jl03308/NAM_pancentromere/genome_alignment/chr8_shujun'
# gff_file1 <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.gff'
# gff_file2 <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/P39.pseudomolecules-v2.fasta.mod.EDTA.TEanno.gff'
# centromere_file <- '/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum'
# genome_size_file <- 'B73.PLATINUM.pseudomolecules-v1.chrom.sizes'
# paf_file <- '/scratch/jl03308/NAM_pancentromere/genome_alignment/chr8_shujun/P39_chr8.mapped-to.B73_chr8.sorted.noseq.paf'
# #core_file <- 'B73_NAM.core.bed'
# chr <- 'chr8'
#alignment_file <- '/scratch/jl03308/NAM_pancentromere/NAM_SV/chr8/B73_P39.aligned.bed'
# gap_file <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/NAM.100ngap.bed'
# ref<-'B73'
# query<-'P39'

teosinte_wdir<-'/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/'
reffiles <- list.files(teosinte_wdir,pattern = paste(ref,".RPKM.bedgraph$",sep=""), include.dirs = TRUE,recursive = TRUE,full.names = TRUE)
queryfiles <- list.files(teosinte_wdir,pattern = paste(query,".RPKM.bedgraph$",sep=""), include.dirs = TRUE,recursive = TRUE,full.names = TRUE)
reffiles_parv<- reffiles[-c(1,2,11,20)]
reffiles_mexicana<- reffiles[c(1,2,11,20)]
queryfiles_parv<- queryfiles[-c(1,2,11,20)]
queryfiles_mexicana<- queryfiles[c(1,2,11,20)]

options(bitmapType='cairo')

library(karyoploteR)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
setwd(wd_dir)
pp <- getDefaultPlotParams(plot.type =1)
pp$data1height <- 400

pp$ideogramheight <- 0.5
pp$data1inmargin<-10
pp$data2inmargin<-10

#gff_file<-'/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.gff'
#chr<-'chr2'
gff<-read.table(gff_file1)
ref.knob180 <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="knob/knob180"),c('V1','V4','V5')]
ref.TR1 <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="knob/TR-1"),c('V1','V4','V5')]
ref.subtelomere <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="subtelomere/4-12-1"),c('V1','V4','V5')]
ref.CentC <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="Cent/CentC"),c('V1','V4','V5')]

if(nrow(ref.TR1) > 0){ref.TR1$V1 <-ref}
if(nrow(ref.knob180) > 0){ref.knob180$V1 <-ref}
if(nrow(ref.subtelomere) > 0){ref.subtelomere$V1 <-ref}
if(nrow(ref.CentC) > 0){ref.CentC$V1 <-ref}
gff<-read.table(gff_file2)
query.knob180 <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="knob/knob180"),c('V1','V4','V5')]
query.TR1 <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="knob/TR-1"),c('V1','V4','V5')]
query.subtelomere <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="subtelomere/4-12-1"),c('V1','V4','V5')]
query.CentC <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="Cent/CentC"),c('V1','V4','V5')]
if(nrow(query.TR1) > 0){query.TR1$V1 <-query}
if(nrow(query.knob180) > 0){query.knob180$V1 <-query}
if(nrow(query.subtelomere) > 0){query.subtelomere$V1 <-query}
if(nrow(query.CentC) > 0){query.CentC$V1 <-query}

#centromere_file<-'/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum'
centromere <- read.table(centromere_file)
ref.centromere <- centromere[(centromere$V1==ref)&(centromere$V2==chr),c('V1','V3','V4')]
query.centromere <- centromere[(centromere$V1==query)&(centromere$V2==chr),c('V1','V3','V4')]

ngap <- read.table(gap_file)
ngap <- ngap[(ngap$V2==chr),c('V1','V3','V4')]


paf<-read.table(paf_file)
genome_size1<-unique(paf$V7)
genome_size2<-unique(paf$V2)

genome <- data.frame("line" = c(ref), "Start" = c(0), "End" = c(genome_size1))
png(paste(ref,".",chr,".teosinte.pairwise.syn.png",sep=""),width=12, height=10, units="in", res=500)
#pdf(paste(ref,"_",query,".",chr,".40-60.pairwise.syn.pdf",sep=""), width=10, height=3)

chromosome <- plotKaryotype(genome = genome, plot.type=1, plot.params = pp)

#palette <- c("#e7eaed","orange","purple","black")
start=0
end=0.05
for(x in reffiles_mexicana){
  print(x)
  #x<-'/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/TIL06/TIL06.B73.RPKM.bedgraph'
  bedgraph<-read.table(x)
  bedgraph_chr=bedgraph[(bedgraph$V1==chr),]
  bedgraph_chr$V1<-ref
  bedgraph_chr[bedgraph_chr$V4>1.5,]$V4 <- 1.5
  kpHeatmap(chromosome, toGRanges(bedgraph_chr), y=bedgraph_chr$V4, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#6d1c49',"#522b6c","black"), ymin=0, ymax=1.5, r0=start,r1=end)
  prefix<-strsplit(basename(x),"[.]")[[1]][1]
  #kpRect(chromosome, toGRanges(line_specific_M1), x0=line_specific_M1$V3,x1=line_specific_M1$V4,y0=0, y1=1,lwd=0.2,border="#6d81c7",color="#6d81c7", r0=start,r1=end)
  #kpPoints(chromosome, chr='B73', x=line_specific$V3, y=line_specific$V5, ymin=0, ymax=max(div[!is.na(div$V5),]$V5), r0=start,r1=end)
  kpAddLabels(chromosome, labels=prefix, r0=start, r1=end,side="left",label.margin=0.005,cex=.9)
  start=start+0.05
  end=end+0.05
}

for(x in reffiles_parv){
  print(x)
  #x<-'/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/TIL06/TIL06.B73.RPKM.bedgraph'
  bedgraph<-read.table(x)
  bedgraph_chr=bedgraph[(bedgraph$V1==chr),]
  bedgraph_chr$V1<-ref
  bedgraph_chr[bedgraph_chr$V4>1.5,]$V4 <- 1.5
  kpHeatmap(chromosome, toGRanges(bedgraph_chr), y=bedgraph_chr$V4, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#6d1c49',"#522b6c","black"), ymin=0, ymax=1.5, r0=start,r1=end)
  prefix<-strsplit(basename(x),"[.]")[[1]][1]
  #kpRect(chromosome, toGRanges(line_specific_M1), x0=line_specific_M1$V3,x1=line_specific_M1$V4,y0=0, y1=1,lwd=0.2,border="#6d81c7",color="#6d81c7", r0=start,r1=end)
  #kpPoints(chromosome, chr='B73', x=line_specific$V3, y=line_specific$V5, ymin=0, ymax=max(div[!is.na(div$V5),]$V5), r0=start,r1=end)
  kpAddLabels(chromosome, labels=prefix, r0=start, r1=end,side="left",label.margin=0.005,cex=.9)
  start=start+0.05
  end=end+0.05
}

kpSegments(chromosome, data=toGRanges(ngap),x0=ngap$V3,x1=ngap$V4,y0=0, y1=1, col="black", 
           lty=5, r0=0, r1=-0.2)

kpRect(chromosome, data=toGRanges(ref.centromere),x0=ref.centromere$V3,x1=ref.centromere$V4,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='#ECCBAE',border='black')
kpRect(chromosome, data=toGRanges(ref.CentC),x0=ref.CentC$V4,x1=ref.CentC$V5,y0=0, y1=1,r0=0, r1=-0.15,col='orange',border=NA)
kpRect(chromosome, data=toGRanges(ref.knob180),x0=ref.knob180$V4,x1=ref.knob180$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='#0000CD',border=NA)
kpRect(chromosome, data=toGRanges(ref.TR1),x0=ref.TR1$V4,x1=ref.TR1$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='#B22222',border=NA)
kpRect(chromosome, data=toGRanges(ref.subtelomere),x0=ref.subtelomere$V4,x1=ref.subtelomere$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='black',border=NA)

kpAddBaseNumbers(chromosome, tick.dist = 10000000, add.units = FALSE,digits=1,cex=0,tick.len = 15,minor.tick.dist = 1000000, minor.tick.len = 5,clipping=TRUE)

dev.off()


genome <- data.frame("line" = c(query), "Start" = c(0), "End" = c(genome_size1))
png(paste(query,".",chr,".teosinte.pairwise.syn.png",sep=""),width=12, height=10, units="in", res=500)

chromosome <- plotKaryotype(genome = genome, plot.type=1, plot.params = pp)

#palette <- c("#e7eaed","orange","purple","black")
start=0
end=0.05
for(x in queryfiles_mexicana){
  print(x)
  #x<-'/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/TIL06/TIL06.B73.RPKM.bedgraph'
  bedgraph<-read.table(x)
  bedgraph_chr=bedgraph[(bedgraph$V1==chr),]
  bedgraph_chr$V1<-query
  bedgraph_chr[bedgraph_chr$V4>1.5,]$V4 <- 1.5
  kpHeatmap(chromosome, toGRanges(bedgraph_chr), y=bedgraph_chr$V4, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#6d1c49',"#522b6c","black"), ymin=0, ymax=1.5, r0=start,r1=end)
  prefix<-strsplit(basename(x),"[.]")[[1]][1]
  #kpRect(chromosome, toGRanges(line_specific_M1), x0=line_specific_M1$V3,x1=line_specific_M1$V4,y0=0, y1=1,lwd=0.2,border="#6d81c7",color="#6d81c7", r0=start,r1=end)
  #kpPoints(chromosome, chr='B73', x=line_specific$V3, y=line_specific$V5, ymin=0, ymax=max(div[!is.na(div$V5),]$V5), r0=start,r1=end)
  kpAddLabels(chromosome, labels=prefix, r0=start, r1=end,side="left",label.margin=0.005,cex=.9)
  start=start+0.05
  end=end+0.05
}

for(x in queryfiles_parv){
  print(x)
  #x<-'/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/TIL06/TIL06.B73.RPKM.bedgraph'
  bedgraph<-read.table(x)
  bedgraph_chr=bedgraph[(bedgraph$V1==chr),]
  bedgraph_chr$V1<-query
  bedgraph_chr[bedgraph_chr$V4>1.5,]$V4 <- 1.5
  kpHeatmap(chromosome, toGRanges(bedgraph_chr), y=bedgraph_chr$V4, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#6d1c49',"#522b6c","black"), ymin=0, ymax=1.5, r0=start,r1=end)
  prefix<-strsplit(basename(x),"[.]")[[1]][1]
  #kpRect(chromosome, toGRanges(line_specific_M1), x0=line_specific_M1$V3,x1=line_specific_M1$V4,y0=0, y1=1,lwd=0.2,border="#6d81c7",color="#6d81c7", r0=start,r1=end)
  #kpPoints(chromosome, chr='B73', x=line_specific$V3, y=line_specific$V5, ymin=0, ymax=max(div[!is.na(div$V5),]$V5), r0=start,r1=end)
  kpAddLabels(chromosome, labels=prefix, r0=start, r1=end,side="left",label.margin=0.005,cex=.9)
  start=start+0.05
  end=end+0.05
}

kpSegments(chromosome, data=toGRanges(ngap),x0=ngap$V3,x1=ngap$V4,y0=0, y1=1, col="black", 
           lty=5, r0=0, r1=-0.2)

kpRect(chromosome, data=toGRanges(query.centromere),x0=query.centromere$V3,x1=query.centromere$V4,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='#ECCBAE',border='black')
kpRect(chromosome, data=toGRanges(query.CentC),x0=query.CentC$V4,x1=query.CentC$V5,y0=0, y1=1,r0=0, r1=-0.15,col='orange',border=NA)
kpRect(chromosome, data=toGRanges(query.knob180),x0=query.knob180$V4,x1=query.knob180$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='#0000CD',border=NA)
kpRect(chromosome, data=toGRanges(query.TR1),x0=query.TR1$V4,x1=query.TR1$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='#B22222',border=NA)
kpRect(chromosome, data=toGRanges(query.subtelomere),x0=query.subtelomere$V4,x1=query.subtelomere$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='black',border=NA)

kpAddBaseNumbers(chromosome, tick.dist = 10000000, add.units = FALSE,digits=1,cex=0,tick.len = 15,minor.tick.dist = 1000000, minor.tick.len = 5,clipping=TRUE)

dev.off()




genome <- data.frame("line" = c(ref,query), "Start" = c(0,0), "End" = c(genome_size1,genome_size2))


alignment_extract <- function(alignment,ref,query) {
#query<-unique(alignment$V4)
alignment1<-alignment[alignment$V10 =="syn",]
#alignment1<-alignment1[(alignment1$V2>40000000)&(alignment1$V3<60000000)&(alignment1$V5>40000000)&(alignment1$V6<60000000),]
data_r<-alignment1[,c('V1','V2','V3','V7')]
data_q<-alignment1[,c('V4','V5','V6')]
data_r_equal<-data_r[data_r$V2==data_r$V3,]
data_r_reverse<-data_r[data_r$V2>data_r$V3,]
data_r_forward<-data_r[data_r$V2<data_r$V3,]
data_r_forward2<-data_r_reverse[,c('V1','V3','V2','V7')]
data_r_equal$V3<-data_r_equal$V2+1
data_r<-rbind(data_r_reverse,data_r_forward2,data_r_equal)
data_q_equal<-data_q[data_q$V5==data_q$V6,]
data_q_reverse<-data_q[data_q$V5>data_q$V6,]
data_q_forward<-data_q[data_q$V5<data_q$V6,]
data_q_forward2<-data_q_reverse[,c('V4','V6','V5')]
colnames(data_q_forward2)<-c('V4','V5','V6')
data_q_equal$V6<-data_q_equal$V5+1

data_q<-rbind(data_q_forward,data_q_forward2,data_q_equal)
data_r<-rbind(data_r_forward,data_r_forward2,data_r_equal)

data_r$V1 <-ref
data_q$V4 <-query
data_r1 <-data_r
data_q1 <-data_q

alignment2<-alignment[(alignment$V10 =="inv")&(alignment$V12 > 20000),]
data_r<-alignment2[,c('V1','V2','V3','V7')]
data_q<-alignment2[,c('V4','V5','V6')]
data_r_equal<-data_r[data_r$V2==data_r$V3,]
data_r_reverse<-data_r[data_r$V2>data_r$V3,]
data_r_forward<-data_r[data_r$V2<data_r$V3,]
data_r_forward2<-data_r_reverse[,c('V1','V3','V2','V7')]
data_r_equal$V3<-data_r_equal$V2+1
data_r<-rbind(data_r_reverse,data_r_forward2,data_r_equal)
data_q_equal<-data_q[data_q$V5==data_q$V6,]
data_q_reverse<-data_q[data_q$V5>data_q$V6,]
data_q_forward<-data_q[data_q$V5<data_q$V6,]
data_q_forward2<-data_q_reverse[,c('V4','V6','V5')]
colnames(data_q_forward2)<-c('V4','V5','V6')
data_q_equal$V6<-data_q_equal$V5+1

data_q<-rbind(data_q_forward,data_q_forward2,data_q_equal)
data_r<-rbind(data_r_forward,data_r_forward2,data_r_equal)

if(nrow(data_r) > 0){data_r$V1 <-ref}
if(nrow(data_q) > 0){data_q$V4 <-query}
link <- list()
link$fr <- data_r1
link$fq <- data_q1
link$rr <- data_r
link$rq <- data_q

return(link)
}

alignment<-read.table(alignment_file)
links1 <- alignment_extract(alignment,ref,query)

pp <- getDefaultPlotParams(plot.type =2)
pp$data1height <- 400
pp$data2height <- 400
pp$ideogramheight <- 0.5
pp$data1inmargin<-10
pp$data2inmargin<-10

pdf(paste(ref,"_",query,".",chr,".mexicana.pairwise.syn.pdf",sep=""), width=10, height=3)
#pdf(paste(ref,"_",query,".",chr,".40-60.pairwise.syn.pdf",sep=""), width=10, height=3)

chromosome <- plotKaryotype(genome = genome, plot.type=1, plot.params = pp)
kpPlotLinks(chromosome, data=toGRanges(links1$fr), data2=toGRanges(links1$fq),col="#cdd0d5",border=NA) ## original is #A79B8E
kpPlotLinks(chromosome, data=toGRanges(links1$rr), data2=toGRanges(links1$rq),col="#ec726c",border=NA) ## original is #A79B8E

#palette <- c("#e7eaed","orange","purple","black")
start=0.15
end=0.3
for(x in reffiles_mexicana){
  print(x)
  #x<-'/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/TIL06/TIL06.B73.RPKM.bedgraph'
  bedgraph<-read.table(x)
  bedgraph_chr=bedgraph[(bedgraph$V1==chr),]
  bedgraph_chr$V1<-ref
  bedgraph_chr[bedgraph_chr$V4>1.5,]$V4 <- 1.5
  kpHeatmap(chromosome, toGRanges(bedgraph_chr), y=bedgraph_chr$V4, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#6d1c49',"#522b6c","black"), ymin=0, ymax=1.5, r0=start,r1=end)
  prefix<-strsplit(basename(x),"[.]")[[1]][1]
  #kpRect(chromosome, toGRanges(line_specific_M1), x0=line_specific_M1$V3,x1=line_specific_M1$V4,y0=0, y1=1,lwd=0.2,border="#6d81c7",color="#6d81c7", r0=start,r1=end)
  #kpPoints(chromosome, chr='B73', x=line_specific$V3, y=line_specific$V5, ymin=0, ymax=max(div[!is.na(div$V5),]$V5), r0=start,r1=end)
  kpAddLabels(chromosome, labels=prefix, r0=start, r1=end,side="left",label.margin=0.005,cex=.9)
  start=start+0.15
  end=end+0.15
}

start=-0.15
end=-0.3
for(x in queryfiles_mexicana){
  print(x)
  #x<-'/scratch/jl03308/NAM_pancentromere/teosinte/parviglumis/read_alignment/TIL06/TIL06.B73.RPKM.bedgraph'
  bedgraph<-read.table(x)
  bedgraph_chr=bedgraph[(bedgraph$V1==chr),]
  bedgraph_chr$V1<-query
  bedgraph_chr[bedgraph_chr$V4>1.5,]$V4 <- 1.5
  kpHeatmap(chromosome, toGRanges(bedgraph_chr), y=bedgraph_chr$V4, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#6d1c49',"#522b6c","black"), ymin=0, ymax=1.5, r0=start,r1=end)
  prefix<-strsplit(basename(x),"[.]")[[1]][1]
  #kpRect(chromosome, toGRanges(line_specific_M1), x0=line_specific_M1$V3,x1=line_specific_M1$V4,y0=0, y1=1,lwd=0.2,border="#6d81c7",color="#6d81c7", r0=start,r1=end)
  #kpPoints(chromosome, chr='B73', x=line_specific$V3, y=line_specific$V5, ymin=0, ymax=max(div[!is.na(div$V5),]$V5), r0=start,r1=end)
  kpAddLabels(chromosome, labels=prefix, r0=start, r1=end,side="left",label.margin=0.005,cex=.9)
  start=start-0.15
  end=end-0.15
}

kpSegments(chromosome, data=toGRanges(ngap),x0=ngap$V3,x1=ngap$V4,y0=0, y1=1, col="black",
           lty=5, r0=0.3, r1=-0.3)

kpRect(chromosome, data=toGRanges(ref.centromere),x0=ref.centromere$V3,x1=ref.centromere$V4,y0=0, y1=1,lwd=0.00,r0=0, r1=0.15,col='#ECCBAE',border='black')
kpRect(chromosome, data=toGRanges(ref.CentC),x0=ref.CentC$V4,x1=ref.CentC$V5,y0=0, y1=1,r0=0, r1=0.15,col='orange',border=NA)
kpRect(chromosome, data=toGRanges(ref.knob180),x0=ref.knob180$V4,x1=ref.knob180$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.15,col='#0000CD',border=NA)
kpRect(chromosome, data=toGRanges(ref.TR1),x0=ref.TR1$V4,x1=ref.TR1$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.15,col='#B22222',border=NA)
kpRect(chromosome, data=toGRanges(ref.subtelomere),x0=ref.subtelomere$V4,x1=ref.subtelomere$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.15,col='black',border=NA)

kpRect(chromosome, data=toGRanges(query.centromere),x0=query.centromere$V3,x1=query.centromere$V4,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='#ECCBAE',border='black')
kpRect(chromosome, data=toGRanges(query.CentC),x0=query.CentC$V4,x1=query.CentC$V5,y0=0, y1=1,r0=0, r1=-0.15,col='orange',border=NA)
kpRect(chromosome, data=toGRanges(query.knob180),x0=query.knob180$V4,x1=query.knob180$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='#0000CD',border=NA)
kpRect(chromosome, data=toGRanges(query.TR1),x0=query.TR1$V4,x1=query.TR1$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='#B22222',border=NA)
kpRect(chromosome, data=toGRanges(query.subtelomere),x0=query.subtelomere$V4,x1=query.subtelomere$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=-0.15,col='black',border=NA)

kpAddBaseNumbers(chromosome, tick.dist = 10000000, add.units = FALSE,digits=1,cex=0,tick.len = 15,minor.tick.dist = 1000000, minor.tick.len = 5,clipping=TRUE)

dev.off()
