#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
wd_dir <- args[1]
gff_file <- args[2]
centromere_file <- args[3]
genome_size_file <- args[4]
div_file <- args[5]
ref<-args[6]
chr <- args[7]
prefix <- args[8]
start_coord<- args[9]
end_coord<- args[10]

start_coord<- as.numeric(start_coord)
end_coord<- as.numeric(end_coord)

#start_coord<- args[7]
#end_coord<- args[8]

library(karyoploteR)
library(data.table)
library(ggplot2)
library(dplyr) 
library(stringr)
library(tidyr)
options(bitmapType='cairo')

setwd(wd_dir)
pp <- getDefaultPlotParams(plot.type = 2)
pp$ideogramheight <- 0
pp$data1height <- 100
pp$data2height <- 200
pp$data1inmargin<-0
pp$data2inmargin<-5
pp$bottommargin<-0
pp$topmargin<-0
pp$data1outmargin<-5
pp$data2outmargin<-5
pp$rightmargin <- 0.05
pp$leftmargin <- 0.07

#gff_file<-'/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.gff'
#chr<-'chr10'
gff<-read.table(gff_file)
B73.knob180 <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="knob/knob180"),c('V1','V4','V5')]
B73.TR1 <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="knob/TR-1"),c('V1','V4','V5')]
B73.subtelomere <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="subtelomere/4-12-1"),c('V1','V4','V5')]
B73.CentC <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="Cent/CentC"),c('V1','V4','V5')]
B73.nor <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="rDNA/spacer"),c('V1','V4','V5')]

#centromere_file<-'/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum'
centromere <- read.table(centromere_file)
B73.centromere <- centromere[(centromere$V1==ref)&(centromere$V2==chr),c('V1','V3','V4')]
#B73.TR1.chr<-B73.TR1[B73.TR1$V1==chr,]
if(nrow(B73.TR1) > 0){B73.TR1$V1 <-ref}
#B73.knob180.chr<-B73.knob180[B73.knob180$V1==chr,]
if(nrow(B73.knob180) > 0){B73.knob180$V1 <-ref}
#B73.subtelomere.chr<-B73.subtelomere[B73.subtelomere$V1==chr,]
if(nrow(B73.subtelomere) > 0){B73.subtelomere$V1 <-ref}
#B73.CentC.chr<-B73.CentC[B73.CentC$V1==chr,]
if(nrow(B73.CentC) > 0){B73.CentC$V1 <-ref}
if(nrow(B73.nor) > 0){B73.nor$V1 <-ref}

#wd_dir<-'/Users/jianingliu/Downloads'
#div_file<-'B73.chr3.div.aligned.txt'
#div_file<-'/scratch/jl03308/NAM_pancentromere/NAM_SV/teosinte/Teo_NAM.P39.flt.vcf.bed'

div<-read.table(div_file)
#tail(div)
query=unique(div[div$V6=="NAM",]$V1)
query2=unique(div[div$V6=="Teo",]$V1)

#genome_size_file<-'/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/P39.pseudomolecules-v2.chrom.sizes'
genome_size<-read.table(genome_size_file)
genome<-data.frame("line" = c(ref), "Start" = c(start_coord), "End" = c(end_coord))
#genome<-data.frame("line" = c(ref), "Start" = c(0), "End" = c(genome_size[genome_size$V1==chr,]$V2))
#genome<-data.frame("line" = c('B73'), "Start" = c(0), "End" = c(182411202))
#genome<-data.frame("line" = c('B73'), "Start" = c(0), "End" = c(243673018))

png(paste(chr,prefix,start_coord,end_coord,"normalized.div.png",sep="."), width=13, height=10, units="in", res=500)
#pdf(paste(chr,".div.pdf",sep=""), width = 10, height = 6)

chromosome <- plotKaryotype(genome = genome, plot.type=2, plot.params = pp)
start=-0.15
end=-0.2
for(x in query){
  #print(x)
  #x<-'B73'
  line_specific=div[(div$V1==x)&(div$V2==chr)&(div$V3>=start_coord)&(div$V3<=end_coord),c('V1','V3','V4','V5')]
  line_specific$V1<-ref
  if(nrow(line_specific[line_specific$V5>200,]) > 0){line_specific[line_specific$V5>200,]$V5 = 200}
  #line_specific=head(line_specific)
  #quantile(div[!is.na(div$V5),]$V5, probs = 0.02,names = FALSE)
  kpHeatmap(chromosome, toGRanges(line_specific), y=line_specific$V5, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#cc0000','#6d1c49',"#522b6c","black"), ymin=0, ymax=200, r0=start,r1=end,clipping=TRUE)
  #kpHeatmap(chromosome, toGRanges(line_specific), y=line_specific$V7, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#6d1c49',"#522b6c","black"), ymin=min(div[!is.na(div$V7),]$V7), ymax=max(max(div[!is.na(div$V7),]$V7),500000), r0=start,r1=end)
  #kpPoints(chromosome, chr='B73', x=line_specific$V3, y=line_specific$V5, ymin=0, ymax=max(div[!is.na(div$V5),]$V5), r0=start,r1=end)
  kpAddLabels(chromosome, labels=x, r0=start, r1=end,side="left",label.margin=0.002,cex=.9)
  start=start-0.05
  end=end-0.05
}

start=0.05
end=0.1
for(x in query2){
  #print(x)
  #x<-'B73'
  line_specific=div[(div$V1==x)&(div$V2==chr)&(div$V3>=start_coord)&(div$V3<=end_coord),c('V1','V3','V4','V5')]
  line_specific$V1<-ref
  if(nrow(line_specific[line_specific$V5>200,]) > 0){line_specific[line_specific$V5>200,]$V5 = 200}
  #line_specific=head(line_specific)
  #quantile(div[!is.na(div$V5),]$V5, probs = 0.02,names = FALSE)
  kpHeatmap(chromosome, toGRanges(line_specific), y=line_specific$V5, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#cc0000','#6d1c49',"#522b6c","black"), ymin=0, ymax=200, r0=start,r1=end,clipping=TRUE)
  #kpHeatmap(chromosome, toGRanges(line_specific), y=line_specific$V7, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#6d1c49',"#522b6c","black"), ymin=min(div[!is.na(div$V7),]$V7), ymax=max(max(div[!is.na(div$V7),]$V7),500000), r0=start,r1=end)
  #kpPoints(chromosome, chr='B73', x=line_specific$V3, y=line_specific$V5, ymin=0, ymax=max(div[!is.na(div$V5),]$V5), r0=start,r1=end)
  kpAddLabels(chromosome, labels=x, r0=start, r1=end,side="left",label.margin=0.002,cex=.9)
  start=start+0.05
  end=end+0.05
}

#head(B73.CentC)
#rownames(B73.CentC) <- NULL
#head(B73.knob180)
#if(nrow(core)>0){kpRect(chromosome, data=toGRanges(core),x0=core$V2,x1=core$V3,y0=0, y1=1,lwd=0.00,r0=0.05, r1=0.1,col='black',border=NA)}
kpRect(chromosome, data=toGRanges(B73.centromere),x0=B73.centromere$V3,x1=B73.centromere$V4,y0=0, y1=1,lwd=0.00,r0=0, r1=0.05,col='#ECCBAE',border='black')
kpRect(chromosome, data=toGRanges(B73.CentC),x0=B73.CentC$V4,x1=B73.CentC$V5,y0=0, y1=1,r0=0, r1=0.05,col='orange',border=NA)
kpRect(chromosome, data=toGRanges(B73.knob180),x0=B73.knob180$V4,x1=B73.knob180$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.05,col='#0000CD',border=NA)
kpRect(chromosome, data=toGRanges(B73.TR1),x0=B73.TR1$V4,x1=B73.TR1$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.05,col='#B22222',border=NA)
kpRect(chromosome, data=toGRanges(B73.subtelomere),x0=B73.subtelomere$V4,x1=B73.subtelomere$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.05,col='black',border=NA)
kpRect(chromosome, data=toGRanges(B73.nor),x0=B73.nor$V4,x1=B73.nor$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.05,col='#198a87',border=NA)
kpAddBaseNumbers(chromosome, tick.dist = 20000000, add.units = FALSE,digits=1,cex=0,tick.len = 6,minor.tick.dist = 1000000, minor.tick.len = 3)

dev.off()


