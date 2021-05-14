library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(reshape2)
#library(fuzzySim)
#require(FactoMineR)
library(cowplot)


options(bitmapType='cairo')


saturation_stat_file='/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/B73.alignment_saturation.bed'
stat <- data.frame(read.table(saturation_stat_file))
genome_file<-'/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/NAM.genome'
genome <- data.frame(read.table(genome_file))
genome<-genome[genome$V2=='B73',]

stat_reshape <- melt(stat,id=c('V1','V2'))

chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')

for(selected_chr in chr){
  stat_reshape$size[stat_reshape$V1==selected_chr] <- genome[genome$V1==selected_chr,]$V4
}

stat_reshape$fraction <- stat_reshape$value/stat_reshape$size 

stat_reshape$V1 <- factor(stat_reshape$V1,
                            levels = chr)


pdf('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/saturation_B73.pdf', width=12, height=15)
ggplot(stat_reshape, aes(x = V2, y = fraction,group=V2)) +
  geom_boxplot() +
  stat_summary(fun=mean, geom="line", aes(group=1))+ 
  stat_summary(fun=mean, geom="point") + 
  xlab("Query") + 
  ylab("Genome aligned with B73 (%)") + 
  scale_x_continuous(breaks = unique(stat_reshape$V2)) +
  #scale_y_continuous(limits = c(0,end)) +
  #scale_y_continuous(trans='log10') + 
  #scale_color_manual(values=col_vector) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=20,angle = 40, hjust = 1),
        axis.text.y=element_text(size=15),
        axis.title=element_text(size=24,face="bold"), 
        strip.text.x = element_text(size =16, colour = "black"),
        strip.text.y = element_text(size =20, colour = "black"),
        legend.text=element_text(size=22),
        legend.title=element_text(size=22)) +
  facet_grid(V1~.) 
dev.off()


stat_reshape_sum <-aggregate(stat_reshape$value, by=list(V2=stat_reshape$V2,variable=stat_reshape$variable), FUN=sum)
#stat_reshape$fraction <- stat_reshape$x/sum(genome$V4)

# Can save the plot with ggsave()
pdf('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/saturation_B73_sum.pdf', width=15, height=12)
ggplot(stat_reshape_sum, aes(x = V2, y = x/sum(genome$V4),group=V2)) +
  geom_boxplot() +
  stat_summary(fun=mean, geom="line", aes(group=1))+ 
  stat_summary(fun=mean, geom="point") + 
  xlab("Query") + 
  ylab("Genome aligned with B73 (%)") + 
  scale_x_continuous(breaks = unique(stat_reshape_sum$V2)) +
  scale_y_continuous(labels=function(x)x,sec.axis = sec_axis(~ . *sum(genome$V4)/1000000000, name = "Size (Gb)")) + 
  #scale_y_continuous(limits = c(0,end)) +
  #scale_y_continuous(trans='log10') + 
  #scale_color_manual(values=col_vector) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=24),
        axis.text.y=element_text(size=24),
        axis.title=element_text(size=28,face="bold"), 
        strip.text.x = element_text(size =20, colour = "black"),
        strip.text.y = element_text(size =24, colour = "black"),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24))
dev.off()
