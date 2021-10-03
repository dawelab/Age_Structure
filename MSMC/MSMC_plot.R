#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
msmc2_file <- args[1]

options(bitmapType='cairo')
library(ggplot2)
require(dplyr)   # for data construction
require(scales)  # for modifying the y-axis

mu <- 3.3e-8
gen <- 1

#msmc2_file<-"/scratch/jl03308/NAM_pancentromere/teosinte_divergence/MSMC/NAM/chr10/NAM.chr10.msmc2.final.txt"
#prefix<-'/scratch/jl03308/NAM_pancentromere/teosinte_divergence/MSMC/NAM/chr10/NAM.chr10.msmc2'

prefix<-file.path(dirname(msmc2_file),paste(unlist(strsplit(basename(msmc2_file),'\\.'))[1:3],collapse="."))
NAM<-read.table(msmc2_file, header=TRUE)
NAM_plot<-data.frame()
Div<-NAM$left_time_boundary/mu*gen
popSize<- (1/NAM$lambda)/mu
NAM_plot <- data.frame(Div,popSize)
breaks <- 10^(1:6)
minor_breaks <- rep(1:9, 21)*(10^rep(1:6, each=9))

png(paste(prefix,"png",sep="."), width=8, height=7, units="in", res=500)
ggplot(data = NAM_plot, mapping = aes(x = Div, y = popSize)) +
  geom_step()+
  scale_y_continuous(trans='log10')+
  scale_x_log10(limits = c(1e1,1e6),breaks = breaks, minor_breaks = minor_breaks)+
  ylab("Effective population size") +
  xlab("Years ago") +
  labs(fill = "") +
  #scale_x_continuous(breaks = unique(stat_umr_plot$index)) +
  #scale_fill_manual(values=c("#334863", "#eda64a")) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=20,angle = 40, hjust = 1),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=24,face="bold"),
        strip.text.x = element_text(size =16, colour = "black"),
        strip.text.y = element_text(size =20, colour = "black"),
        legend.text=element_text(size=22),
        legend.title=element_text(size=22))

dev.off()
