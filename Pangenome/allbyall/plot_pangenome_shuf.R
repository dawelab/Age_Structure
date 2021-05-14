library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(reshape2)

options(bitmapType='cairo')


saturation_stat_file='/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/pangenome/wholegenome.pangenome.bed'
saturation_stat <- data.frame(read.table(saturation_stat_file))
#saturation_stat_reshape <- melt(saturation_stat[,c('V1','V3')],id=c('V1'))
colnames(saturation_stat) <- c('index','line','size','permutation','chr')
#stat_umr_plot<-melt(stat_umr, id=c("index",'chr'))
chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
saturation_stat$chr <- factor(saturation_stat$chr,
                            levels = chr)

pdf('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/pangenome/pangenome.pdf', width=12, height=15)

ggplot(saturation_stat, aes(x = index, y = size/1000000,group=index)) +
  geom_boxplot() +
  stat_summary(fun=mean, geom="line", aes(group=1))+ 
  stat_summary(fun=mean, geom="point") + 
  xlab("Genomes") + 
  ylab("Pangenome space (Mb)") + 
  scale_x_discrete(limits = unique(saturation_stat$index)) +
  #scale_y_continuous(limits = c(0,end)) +
  #scale_y_continuous(trans='log10') + 
  #scale_color_manual(values=col_vector) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=20,angle = 40, hjust = 1),
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=24,face="bold"), 
        strip.text.x = element_text(size =16, colour = "black"),
        strip.text.y = element_text(size =20, colour = "black"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20)) + 
  facet_grid(chr~.,scales="free_y") 

dev.off()

saturation_stat_sum<-aggregate(saturation_stat$size, by=list(index=saturation_stat$index,permutation=saturation_stat$permutation), FUN=sum)
colnames(saturation_stat_sum) <- c('index','permutation','size')
pangenomesize<-mean(saturation_stat_sum[saturation_stat_sum$index==26,]$size)/1000000000
pdf('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/pangenome/pangenome.wholegenome.pdf', width=12, height=8)

ggplot(saturation_stat_sum, aes(x = index, y = size/1000000000,group=index)) +
  geom_boxplot() +
  stat_summary(fun=mean, geom="line", aes(group=1))+ 
  scale_y_continuous(labels=function(x)x,sec.axis = sec_axis(~ . /pangenomesize, name = "Fraction"))+ 
  stat_summary(fun=mean, geom="point") + 
  xlab("Genomes") + 
  ylab("Pangenome space (Gb)") + 
  scale_x_discrete(limits = unique(saturation_stat_sum$index)) +
  #scale_y_continuous(limits = c(0,end)) +
  #scale_y_continuous(trans='log10') + 
  #scale_color_manual(values=col_vector) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=20,angle = 40, hjust = 1),
        axis.text.y=element_text(size=20),
        axis.title=element_text(size=24,face="bold"), 
        strip.text.x = element_text(size =16, colour = "black"),
        strip.text.y = element_text(size =20, colour = "black"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20)) 

dev.off()
