#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
wdir <- args[1]
snp_file <- args[2]
output_file <- args[3]

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(plyr)
library(reshape2)

options(bitmapType='cairo')

wdir<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/'
snp_file<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/NAM.knob.dating.bed'
output_file<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/NAM.knob.dating'
setwd(wdir)
nam_snp <- data.frame(read.table(snp_file))
nam_snp_plot <- nam_snp[,c('V4','V15','V16')]
#colnames(nam_centromere) <- c("NAM_line", "Chr","Start","End","Size","100N gap number","13N gap number","Known gap size","Fully assembled (Y/n)","CentC","CRM2","CRM1","cinful-zeon","opie-ji","huck","prem1","grande", "All_TE","gene","Known sequences", "Other TE","Total gap","Unknown" )
colnames(nam_snp_plot) <- c("NAM_line", 'Divergence','Knob_type')
div_mean <- ddply(nam_snp, .(V4,V16), summarize, Divergence = mean(V13/V14/2/3.3*100))
colnames(div_mean) <- c("NAM_line", 'Knob_type', 'Divergence')

#nam_centromere_plot$Chr <- factor(nam_centromere_plot$Chr,levels = chr)
png(paste(output_file,"distribution.png",sep="."), width=10, height=12, units="in", res=500)

ggplot(nam_snp_plot, aes(x=NAM_line, y=Divergence/1000000)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(size=0.8,alpha=0.6) + 
  geom_point(data = div_mean, 
             mapping = aes(x = NAM_line, y = Divergence),
             col="red") +
  geom_line(data = div_mean, 
            mapping = aes(x = NAM_line, y = Divergence, group=Knob_type, color='red')) + 
  #stat_summary(fun.y=median, geom="point", size=2) + 
  #geom_boxplot(width=0.1) +
  ylab("Divergence Time (million years)") + 
  xlab("Line") + 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values=c("#87a2ba", "#334863", "#eda64a","#90766b","#ead09e","#ec853f","#00a950","#1c8b82","#ab6661","#79538e","#a2a5b4")) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) +
  geom_text(data = div_mean, aes(y = Divergence, label = round(Divergence,2)),size = 3.5, col="red",vjust=-5,position = position_dodge(width = 1)) +
  scale_y_continuous(limits = c(0,1.5)) + 
  facet_grid(Knob_type~.) 
dev.off()

