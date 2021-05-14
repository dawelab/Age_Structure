#!/usr/bin/env Rscript
options(bitmapType='cairo')

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)

wdir<-'/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis'
allelefrequency_file='/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/alignment_frequency_allchrs.bed'

# plot the size distribution of each frequency 
allelefrequency_plot$size <- allelefrequency_plot$end- allelefrequency_plot$start +1 

png('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency_size.png', width=12, height=18, units="in", res=500)
ggplot(allelefrequency_plot, aes(y=size, x=sum)) + 
  
  geom_jitter(position=position_jitter(0.2), cex=0.2) + 
  #geom_boxplot(alpha=0.8) + 
  scale_y_continuous(trans='log10',labels = scales::comma)+
  stat_summary(fun.y=median, geom = "crossbar",width = .5, color = "red") +
  ylab("Size of fragments (%)") + 
  xlab("Alle frequency among 25 NAM genomes") + 
  labs(fill = "") +
  scale_x_discrete(breaks = sort(unique(allelefrequency_plot$sum))) +
  #scale_fill_manual(values=c("#334863", "#eda64a")) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=20,angle = 40, hjust = 1),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=24,face="bold"), 
        strip.text.x = element_text(size =16, colour = "black"),
        strip.text.y = element_text(size =20, colour = "black"),
        legend.text=element_text(size=22),
        legend.title=element_text(size=22)) + 
  facet_grid(chr~.,scales="free_y") 
dev.off()
