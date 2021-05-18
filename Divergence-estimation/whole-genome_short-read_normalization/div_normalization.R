#!/usr/bin/env Rscript
# args <- commandArgs(trailingOnly = TRUE)
# wdir <- args[1]
# snp_file <- args[2]
# output_file <- args[3]

#library(ggplot2)
library(ggpmisc)
options(bitmapType='cairo')
wdir<-'/scratch/jl03308/NAM_pancentromere/teosinte_divergence/'
snp_file<-'/scratch/jl03308/NAM_pancentromere/teosinte_divergence/Teo_NAM.B73.div.bed'
abssnp_file<-'/scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_20000_div.bed'

setwd(wdir)
nam_snp <- data.frame(read.table(snp_file))
nam_snp_nam<- nam_snp[grepl("Maize", nam_snp$V8),]
nam_snp_plot <- nam_snp_nam[,c('V7','V1','V2','V3','V5','V6')]
colnames(nam_snp_plot) <- c("NAM_line", 'Chr','Start','End',"SNP",'Divergence')

nam_abssnp <- data.frame(read.table(abssnp_file))
nam_abssnp_plot <- nam_abssnp[,c('V4','V6','V2','V3','V7','V8')]
colnames(nam_abssnp_plot) <- c("NAM_line", 'Chr','Start','End',"SNP",'Divergence')

nam_comparison <- nam_abssnp_plot %>% right_join(nam_snp_plot, by=c("NAM_line", 'Chr','Start','End'))
nam_comparison$diff <- nam_comparison$Divergence.x/nam_comparison$Divergence.y
nam_comparison$diff <- abs(nam_comparison$Divergence.x - nam_comparison$Divergence.y)
nam_comparison$diff2 <- abs(nam_comparison$Divergence.x/nam_comparison$Divergence.y)

nam_comparison<-nam_comparison[(nam_comparison$SNP.x >0)&(nam_comparison$SNP.x >0),]
png("snp_normalization.png", width=4, height=4, units="in", res=500)
ggplot(nam_comparison, aes(x=Chr, y=diff)) + 
  geom_violin(trim=TRUE) +
  #geom_boxplot(width=0.1) +
  xlab("Chr") + 
  ylab("Fold change") + 
  stat_summary(fun.y=median, geom="point", size=2, color="red")+ 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values=c("#87a2ba", "#334863", "#eda64a","#90766b","#ead09e","#ec853f","#00a950","#1c8b82","#ab6661","#79538e","#a2a5b4")) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) +
  #scale_y_continuous(limits = c(0,5),minor_breaks = seq(0,5, 0.5)) 
  scale_y_continuous(limits = c(0,50000),minor_breaks = seq(0,50000, 5000)) 

dev.off()

library(tidyr)
samplenum <- 10000
nam_comparison_sampled <- sample_n(nam_comparison[nam_comparison$Divergence.x < 1000000,] %>% drop_na(), samplenum)
nam_comparison_sampled$paired <- seq.int(nrow(nam_comparison_sampled))
nam_comparison_sampled_reshape <- melt(nam_comparison_sampled[c('Divergence.x','Divergence.y','paired')],id=c("paired"))
nam_comparison_sampled_reshape$variable <- gsub('Divergence.x', 'syntenic SNPs',nam_comparison_sampled_reshape$variable)
nam_comparison_sampled_reshape$variable <- gsub('Divergence.y', 'short-reads SNPs',nam_comparison_sampled_reshape$variable)

#nam_comparison_sampled_reshape$color <- 0
#nam_comparison_sampled_reshape[(nam_comparison_sampled_reshape$value > 150000) & (nam_comparison_sampled_reshape$variable == 'short-reads SNPs'),]$color <-1
png("snp_normalization.png", width=4, height=4, units="in", res=500)
nam_comparison_sampled_reshape %>%
  ggplot(aes(variable,value,fill=variable)) +
  geom_boxplot() +
  geom_point()+ 
  geom_line(aes(group=paired),size=0.05,alpha=0.1) +
  theme_minimal() +
  xlab("") + 
  ylab("Divergence time (years)") + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_y_continuous(limits = c(0,1000000),minor_breaks = seq(0,1000000, 50000)) 
dev.off()

png("snp_div_distribution.png", width=4, height=4, units="in", res=500)
nam_comparison_sampled %>%
  ggplot(aes(Divergence.x,diff2)) +
  geom_point(cex = 0.2,alpha=0.5) +
  theme_minimal() +
  ylab("Syn SNPs / Short-read SNPs") + 
  xlab("Divergence time (years)") + 
  theme(legend.position = "none") +
  scale_fill_manual(values=c("#999999", "#E69F00"))+
  scale_x_continuous(limits = c(0,1000000),minor_breaks = seq(0,1000000, 50000)) +
  scale_y_continuous(limits = c(0,10)) 
dev.off()

# tripsacum div estimation 
tripsacum_snp<- nam_snp[grepl("tripsacum", nam_snp$V9),]
# 3mb region
tripsacum_snp_nor <- tripsacum_snp[(tripsacum_snp$V1=="chr6")&(tripsacum_snp$V4>50)&(tripsacum_snp$V3<=20960000)&(tripsacum_snp$V2>=17800000),]
# normalization
sum(tripsacum_snp_nor$V5)/sum(tripsacum_snp_nor$V4)/2/3.3 *100000000 * median(sapply(nam_comparison[nam_comparison$Divergence.x > 500000,]$diff2, as.numeric),na.rm=TRUE)

