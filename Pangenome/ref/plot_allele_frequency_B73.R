#!/usr/bin/env Rscript
options(bitmapType='cairo')

library(karyoploteR)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)

pp <- getDefaultPlotParams(plot.type = 1)
pp$data1height <- 400
pp$ideogramheight <- 0.5
pp$data1inmargin<-0

wdir<-'/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis'
genome_file<-'/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/NAM.genome'
allelefrequency_file='/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/alignment_frequency_allchrs.bed'
centromere_file<-'/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum'
gap_file <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/NAM.100ngap.bed'

setwd(wdir)
centromere <- read.table(centromere_file)
#ngap <- read.table(gap_file)
#line.ngap <- ngap[(ngap$V1==line),c('V2','V3','V4')]
genome <- read.table(genome_file)
genome$V3 <-0
allelefrequency <- data.frame(read.table(allelefrequency_file, header=TRUE))
allelefrequency_plot <- allelefrequency[c('chr','start','end','sum')]
chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
allelefrequency_plot$chr <- factor(allelefrequency_plot$chr,
                            levels = chr)

line<-'B73'
line.centromere <- centromere[(centromere$V1==line),c('V2','V3','V4')]
gff_file <- list.files('/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/', pattern=glob2rx(paste(line,'.*pseudomolecules-v*.fasta.mod.EDTA.TEanno.gff',sep='')),full.names=TRUE)
#gff_file <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.gff'
gff<-read.table(gff_file)
line.knob180 <- gff[(gff$V3=="knob/knob180"),c('V1','V4','V5')]  %>% mutate(V1 = gsub(paste(line,"_",sep=''), "", V1))
line.TR1 <- gff[(gff$V3=="knob/TR-1"),c('V1','V4','V5')] %>% mutate(V1 = gsub(paste(line,"_",sep=''), "", V1))
line.subtelomere <- gff[(gff$V3=="subtelomere/4-12-1"),c('V1','V4','V5')] %>% mutate(V1 = gsub(paste(line,"_",sep=''), "", V1))
line.CentC <- gff[(gff$V3=="Cent/CentC"),c('V1','V4','V5')] %>% mutate(V1 = gsub(paste(line,"_",sep=''), "", V1))
line.NOR <- gff[(gff$V3=="rDNA/spacer"),c('V1','V4','V5')] %>% mutate(V1 = gsub(paste(line,"_",sep=''), "", V1))

genomeline<- genome[(genome$V2==line),c('V1','V3','V4')]
#pdf(paste(line,".sv.pdf",sep=""), width=7, height=3)
png(paste(line,".allelefrequency.png",sep=""), width=8, height=6, units="in", res=500)
#pdf(paste(ref,"_",query,".",chr,".40-60.pairwise.syn.pdf",sep=""), width=10, height=3)
chromosome <- plotKaryotype(genome = genomeline, plot.type=1, plot.params = pp)
#kpPlotRegions(chromosome, data=unalignment, col="#cdd0d5", r0=0.05, r1=1,non.overlapping = FALSE, border=NA)
#kpPlotRegions(chromosome, data=inv, col="#ec726c", r0=0.05, r1=1,non.overlapping = FALSE)
#kpPlotRegions(chromosome, data=tandemdup, col="#13477d", r0=0.05, r1=1)
start=0.25
end=0.75
kpHeatmap(chromosome, toGRanges(allelefrequency_plot), y=allelefrequency_plot$sum, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#cc0000',"#522b6c","black"), ymin=0, ymax=25, r0=start,r1=end)

kpRect(chromosome, data=toGRanges(line.centromere),x0=line.centromere$V3,x1=line.centromere$V4,y0=0, y1=1,lwd=0.00,r0=0, r1=0.25,col='#ECCBAE',border='black')
kpRect(chromosome, data=toGRanges(line.CentC),x0=line.CentC$V4,x1=line.CentC$V5,y0=0, y1=1,r0=0, r1=0.25,col='orange',border=NA)
kpRect(chromosome, data=toGRanges(line.knob180),x0=line.knob180$V4,x1=line.knob180$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.25,col='#0000CD',border=NA)
kpRect(chromosome, data=toGRanges(line.TR1),x0=line.TR1$V4,x1=line.TR1$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.25,col='#B22222',border=NA)
kpRect(chromosome, data=toGRanges(line.subtelomere),x0=line.subtelomere$V4,x1=line.subtelomere$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.25,col='black',border=NA)
kpRect(chromosome, data=toGRanges(line.NOR),x0=line.NOR$V4,x1=line.NOR$V5,y0=0, y1=1,lwd=0.00,r0=0, r1=0.25,col='#198a87',border=NA)
kpAddBaseNumbers(chromosome, tick.dist = 20000000, add.units = FALSE,digits=1,cex=0,tick.len = 15,minor.tick.dist = 1000000, minor.tick.len = 5,clipping=TRUE)
dev.off()

library(circlize)
library(ComplexHeatmap)
#pdf(paste("legends.pdf",sep=""), width=5, height=5)
col_fun = colorRamp2(c(0,5,10,15,20,25), c("#e7eaed",'#e6ab00',"#fb7c24",'#cc0000',"#522b6c","black"))
lgd = Legend(col_fun = col_fun, direction = "horizontal",legend_width = unit(8, "cm"), at = c(0,5,10,15,20,25))
draw(lgd, x = unit(0.5, "npc"), y = unit(0.05, "npc"))
dev.off()
