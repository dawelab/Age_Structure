#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
wdir <- args[1]
plot_stat <- args[2]
output<-args[3]

library(reshape2)
library(ggplot2)
library(tidyverse)
library(data.table)
library(facetscales)
library(tools)
library(stringr)

setwd(wdir)

options(bitmapType='cairo')
#plot_stat<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/pairwise_alignment/TR-1/chr9/knob1/TR-1.chr9.knob1.B97.ref.labelled.stat.plot'
align_stat<-read.table(plot_stat)
align_stat <- data.table(align_stat)

type<-unlist(strsplit(file_path_sans_ext(basename(plot_stat)),"\\."))[[1]]
prefix<-unlist(strsplit(file_path_sans_ext(basename(plot_stat)),"\\."))[[3]]
chr<-unlist(strsplit(file_path_sans_ext(basename(plot_stat)),"\\."))[[2]]
repeat_num_file<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/NAM.repeat.num.txt'
repeat_num <- data.frame(read.table(repeat_num_file))
repeat_num_selected <- repeat_num[(repeat_num$V1==type)&(repeat_num$V2==chr)&(repeat_num$V3==prefix),]
nam_lines <- unique(repeat_num_selected$V4)
align_stat$min<-1
align_stat$max<-0
for(x in unique(str_split_fixed(align_stat$V4, "[.]", 2)[,1])){
  align_stat$max[align_stat$V4==x]<-repeat_num_selected[repeat_num_selected$V4==x,]$V5
}


#yaxis_limit<-unique( align_stat[ , c('V4','V5_min','V5_max')] )

align_stat<-align_stat[align_stat$V6>=98,]
align_stat <- align_stat %>% mutate(V2 = gsub(paste(".",prefix,sep=""), "", V2))
align_stat <- align_stat %>% mutate(V4 = gsub(paste(".",prefix,sep=""), "", V4))
ref<-unique(align_stat$V2)
reflen<-repeat_num[(repeat_num$V1==type)&(repeat_num$V2==chr)&(repeat_num$V3==prefix)&(repeat_num$V4==ref),]$V5
#querylen<-repeat_num_selected[repeat_num_selected$V4!=ref,c('V4','V5')]
  
png(output, width=14, height=10, units="in", res=500)
ggplot()+
  geom_point(data = align_stat, aes(V5, V3, color = V6),size = 0.0001,shape=20)+
  scale_color_gradient2(low = "#e5acb6", high = "#0f0300", mid = "#5b3126", midpoint = 99, limit = c(98,100), na.value="white",space = "Lab", name="JI (%)") +
  theme_minimal()+ 
  geom_blank(data = align_stat,aes(x = min)) +
  geom_blank(data = align_stat,aes(x = max)) + 
  #geom_blank(aes(y = y_max))
  #geom_blank(data=querylen,aes(x = V5)) + 
  #coord_fixed() +
  #scale_x_continuous(limits = c(1, x_limit_upper)) + 
  theme(axis.text.x=element_text(size=13,angle = 30, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=14), 
        axis.title.x = element_blank(),
        strip.text.x = element_text(size =13, colour = "black"),
        strip.text.y = element_text(size =13, colour = "black"),
        plot.title = element_text(size=15,face='bold',hjust = 0.5,vjust = 0.5)) + 
        #strip.placement = "outside") +
  facet_wrap(~V4, ncol = 7, scales = "free_x") + 
  coord_cartesian(ylim = c(1, reflen)) +
  ylab(ref)
  
dev.off()


