#!/usr/bin/env Rscript

library(ggplot2)
options(bitmapType='cairo')

#setwd(wdir)
#nam_centromere_plot$Chr <- factor(nam_centromere_plot$Chr,levels = chr)

snp_file <-'/scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_pangenome_div.bed'
prefix<- '/scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_AF_div'
nam_snp <- data.frame(read.table(snp_file))
nam_snp_plot <- nam_snp[,c('V7','V8','V9')]
#colnames(nam_centromere) <- c("NAM_line", "Chr","Start","End","Size","100N gap number","13N gap number","Known gap size","Fully assembled (Y/n)","CentC","CRM2","CRM1","cinful-zeon","opie-ji","huck","prem1","grande", "All_TE","gene","Known sequences", "Other TE","Total gap","Unknown" )
colnames(nam_snp_plot) <- c('Divergence',"chr","index")
chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
nam_snp_plot$chr <- factor(nam_snp_plot$chr,levels = chr)
#nam_snp_plot<-nam_snp_plot[nam_snp_plot$Divergence>50000,]
png(paste(prefix,"boxplot_distribution.png",sep="."), width=8, height=12, units="in", res=500)
ggplot(nam_snp_plot, aes(x=index, y=Divergence/1000000,group=index)) + 
  #geom_violin() + 
  #geom_jitter(position=position_jitter(0.2),cex=0.01,alpha = 0.2) +
  #stat_summary(fun.y=median, geom="point", size=2) + 
  geom_boxplot() +
  ylab("Divergence Time (million years)") + 
  xlab("Allele frequency") + 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values=c("#87a2ba", "#334863", "#eda64a","#90766b","#ead09e","#ec853f","#00a950","#1c8b82","#ab6661","#79538e","#a2a5b4")) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) +
  scale_y_continuous(limits = c(0,1)) + 
  scale_x_discrete(limits=c(0:25)) +
  facet_grid(chr~.) 
dev.off()

nam_snp_plot$index <- paste("AF", nam_snp_plot$index,sep="_")
library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-4]

density_plot <- ggplot(nam_snp_plot, aes(x=Divergence/1000000, color=index)) + 
  #geom_jitter(position=position_jitter(0.2),cex=0.005) +
  geom_density() + 
  #stat_peaks(colour = "red") +
  #stat_valleys(colour = "blue") + 
  #stat_summary(fun.y=median, geom="point", size=2) + 
  #geom_boxplot(width=0.1) +
  xlab("Divergence Time (million years)") + 
  #xlab("Line") + 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values=col_vector) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) +
  scale_x_continuous(limits = c(0.05,1),minor_breaks = seq(0.05,1, 0.05)) #+
#facet_wrap(~ Chr,scales="free", ncol =1, strip.position="right") 
density_plot_build <- ggplot_build(density_plot)

png(paste(prefix,"density_distribution.png",sep="."), width=8, height=12, units="in", res=500)
density_plot + stat_peaks(
  data = density_plot_build[['data']][[1]], # take a look at this object
  aes(x = x, y = density),
  colour = "red",
  size =1.5
) 
#nam_snp_plot_pericentromere<-nam_snp_plot_pericentromere[nam_snp_plot_pericentromere$Divergence >50000,]
dev.off()

