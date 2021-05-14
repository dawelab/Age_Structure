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


stat_file='/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/B73.allele_frequency.stat'
stat <- data.frame(read.table(stat_file))
colnames(stat) <- c('ref','index','totalsize','totalsize2Chrlen','genesize','gene2total','genesize2totalgene','UMRsize','UMR2total','UMR2totalUMR','geneUMRsize','geneUMR2total','geneUMR2totalgeneUMR','chr')
stat$other <- stat$totalsize- stat$geneUMRsize
stat_umr <- stat[,c('index','geneUMRsize','other','chr')]
colnames(stat_umr) <- c('index','Gene&UMR','Other','chr')
stat_umr_plot<-melt(stat_umr, id=c("index",'chr'))
chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
stat_umr_plot$chr <- factor(stat_umr_plot$chr,
                            levels = chr)
for(selected_chr in chr){
  stat_umr_plot$size[stat_umr_plot$chr==selected_chr] <- sum(stat_umr_plot[stat_umr_plot$chr==selected_chr,]$value)
}
stat_umr_plot$fraction <- stat_umr_plot$value/stat_umr_plot$size 
stat_umr_plot$index <- stat_umr_plot$index + 1 
pdf('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency_barplot.pdf', width=12, height=15)
#png('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency_barplot.png', width=12, height=18, units="in", res=500)
# allelefrequency_barplot <- function(data){
  
ggplot(stat_umr_plot, aes(fill=variable, y=fraction, x=index)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_y_continuous(labels=function(x)x*100)+ 
  ylab("Fraction of B73 genome (%)") + 
  xlab("Alle frequency among 26 genomes") + 
  labs(fill = "") +
  scale_x_continuous(breaks = unique(stat_umr_plot$index)) +
  scale_fill_manual(values=c("#334863", "#eda64a")) +
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

#}
#stat_umr_plot2 <- split(stat_umr_plot, stat_umr_plot$chr)
#p_lst <- lapply(stat_umr_plot2, allelefrequency_barplot)

#library(cowplot)
#p0 <- ggplot() + labs(title="Export Ratio and Real Effective Exchange Rate")
#p1 <- plot_grid(plotlist=p_lst, ncol=1)
#pp <- plot_grid(p0, p1, ncol=1, rel_heights=c(0.1, 1))

dev.off()

stat_function <- stat[,c('index','genesize2totalgene','UMR2totalUMR','geneUMR2totalgeneUMR','chr')]
colnames(stat_function) <- c('index','Gene','UMR','Gene&UMR','chr')
stat_function_plot<-melt(stat_function, id=c("index",'chr'))
stat_function_plot$index <- stat_function_plot$index +1
pdf('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency_pieplot.pdf', width=12, height=15)

#png('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency_pieplot.png', width=12, height=18, units="in", res=500)
ggplot(stat_function_plot, aes(fill=as.integer(index), y=value, x=variable)) + 
  geom_col(colour="black") +
  scale_fill_gradient(name = "Allele frequency", 
                      high = "#334863", low = "#eaecef") +
  coord_polar("y") + 
  scale_x_discrete(unique(stat_umr_plot$variable))  +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =element_blank(), 
        axis.text.y=element_text(size =18, colour = "black"),
        strip.text.y = element_text(size =20, colour = "black"),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.key.size = unit(1.2, "cm")) +
  ylab("") + 
  xlab("") +
  facet_grid(chr~.) 
dev.off()

stat_umr <- stat[,c('index','geneUMRsize','other','chr')]
x<-aggregate(stat_umr$geneUMRsize, by=list(index=stat_umr$index), FUN=sum)
y<-aggregate(stat_umr$other, by=list(index=stat_umr$index), FUN=sum)
stat_umr <- merge(x,y,by=c("index"))
colnames(stat_umr) <- c('index','Gene&UMR','Other')
stat_umr_plot<-melt(stat_umr, id=c("index"))
stat_umr_plot$fraction<- stat_umr_plot$value/sum(stat$totalsize)
stat_umr_plot$index <- stat_umr_plot$index + 1 

barplot<-ggplot(stat_umr_plot, aes(fill=variable, y=fraction, x=index)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_y_continuous(labels=function(x)x*100,sec.axis = sec_axis(~ . *sum(stat$totalsize)/1000000, name = "Size (Mb)"))+ 
  ylab("Fraction of B73 genome (%)") + 
  xlab("Alle frequency among 26 genomes") + 
  labs(fill = "") +
  scale_x_continuous(breaks = unique(stat_umr_plot$index)) +
  scale_fill_manual(values=c("#334863", "#eda64a")) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=24,angle = 40, hjust = 1),
        axis.text.y=element_text(size=24),
        axis.title=element_text(size=28,face="bold"), 
        strip.text.x = element_text(size =20, colour = "black"),
        strip.text.y = element_text(size =24, colour = "black"),
        legend.text=element_text(size=24),
        legend.title=element_text(size=24))
  
stat_function <- stat[,c('index','genesize2totalgene','UMR2totalUMR','geneUMR2totalgeneUMR','chr')]
x<-aggregate(stat_function$genesize2totalgene, by=list(index=stat_function$index), FUN=sum)
y<-aggregate(stat_function$UMR2totalUMR, by=list(index=stat_function$index), FUN=sum)
z<-aggregate(stat_function$geneUMR2totalgeneUMR, by=list(index=stat_function$index), FUN=sum)

stat_function <- merge(merge(x,y,by=c("index")),z,by=c("index"))
colnames(stat_function) <- c('index','Gene','UMR','Gene&UMR')
stat_function_plot<-melt(stat_function, id=c("index"))
stat_function_plot$index <- stat_function_plot$index + 1 

pieplot<-ggplot(stat_function_plot, aes(fill=as.integer(index), y=value, x=variable)) + 
  geom_col(colour="black") +
  scale_fill_gradient(name = "Allele frequency", 
                      high = "#334863", low = "#eaecef") +
  coord_polar("y") + 
  scale_x_discrete(unique(stat_umr_plot$variable))  +
  theme_minimal() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y =element_blank(), 
        axis.text.y=element_text(size =25, colour = "black"),
        strip.text.y = element_text(size =25, colour = "black"),
        legend.text=element_text(size=25),
        legend.title=element_text(size=25),
        legend.key.size = unit(1.2, "cm")) +
  ylab("") + 
  xlab("") 

#barplot + annotation_custom(ggplotGrob(pieplot), xmin = 7, xmax = 17, 
#                       ymin = 15, ymax = 26)


# Can save the plot with ggsave()
pdf('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency_sum_barplot.pdf', width=20, height=15)
#png('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency_sum_barplot.png', width=18, height=12, units="in", res=500)
barplot
dev.off()
pdf('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency_sum_pieplot.pdf', width=12, height=15)
#png('/lustre2/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/breakpoint_analysis/allele_frequency_sum_pieplot.png', width=8, height=12, units="in", res=500)
pieplot
dev.off()
