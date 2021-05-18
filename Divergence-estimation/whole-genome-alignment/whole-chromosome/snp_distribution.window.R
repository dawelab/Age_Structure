#!/usr/bin/env Rscript
# args <- commandArgs(trailingOnly = TRUE)
# wdir <- args[1]
# snp_file <- args[2]
# output_file <- args[3]

#library(ggplot2)
library(ggpmisc)

options(bitmapType='cairo')
wdir<-'/scratch/jl03308/NAM_pancentromere/divergence/intergenic/'
snp_file<-'/scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_20000_div.bed'
prefix<- '/scratch/jl03308/NAM_pancentromere/divergence/intergenic/NAM.intergenic_20000_div'
setwd(wdir)

nam_snp <- data.frame(read.table(snp_file))
nam_snp_plot <- nam_snp[,c('V4','V6','V2','V3','V7','V8')]
#colnames(nam_centromere) <- c("NAM_line", "Chr","Start","End","Size","100N gap number","13N gap number","Known gap size","Fully assembled (Y/n)","CentC","CRM2","CRM1","cinful-zeon","opie-ji","huck","prem1","grande", "All_TE","gene","Known sequences", "Other TE","Total gap","Unknown" )
colnames(nam_snp_plot) <- c("NAM_line", 'Chr','Start','End',"SNP",'Divergence')
chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
nam_snp_plot$Chr <- factor(nam_snp_plot$Chr,
                           levels = chr)

#nam_centromere_plot$Chr <- factor(nam_centromere_plot$Chr,levels = chr)
# png(paste(output_file,"distribution.png",sep="."), width=8, height=15, units="in", res=500)
# 
# ggplot(nam_snp_plot, aes(x=NAM_line, y=Divergence/1000000)) + 
#   geom_violin() + 
#   #stat_summary(fun.y=median, geom="point", size=2) + 
#   #geom_boxplot(width=0.1) +
#   ylab("Divergence Time (million years)") + 
#   xlab("Line") + 
#   #scale_y_continuous(trans='log10') + 
#   scale_color_manual(values=c("#87a2ba", "#334863", "#eda64a","#90766b","#ead09e","#ec853f","#00a950","#1c8b82","#ab6661","#79538e","#a2a5b4")) +
#   #theme(legend.position="top") +
#   theme_minimal() +
#   theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
#         axis.text.y=element_text(size=13),
#         axis.title=element_text(size=16,face="bold"), 
#         strip.text.x = element_text(size =15, colour = "black"),
#         strip.text.y = element_text(size =15, colour = "black")) +
#   scale_y_continuous(limits = c(0,1)) + 
#   facet_grid(Chr~.) 
# dev.off()

library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[-4]

png(paste(prefix,"wholechr.png",sep="."), width=15, height=22, units="in", res=500)
ggplot(nam_snp_plot, aes(x=Start/1000000, y=Divergence/1000000, color=NAM_line)) +
  geom_point(cex = 0.005,alpha=0.7) + 
  ylab("Divergence Time (million years)") + 
  xlab("Chromosome (Mb)") + 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values=col_vector) +
  #theme(legend.position="top") +
  guides(color=guide_legend(override.aes = list(size = 3,alpha=0.8),ncol = 9)) + 
  theme_minimal() +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=18),
        axis.title=element_text(size=24,face="bold"),
        strip.text = element_text(size =24, colour = "black"),
        legend.title = element_blank(),
        legend.position="bottom",legend.text = element_text(size=16)) +
  scale_y_continuous(limits = c(0,1)) + 
  #facet_wrap(~ Chr,scales="free", ncol =1, strip.position="right") 
  facet_grid(Chr~.,scales="free")
dev.off()




png(paste(prefix,"wholechr.density_distribution.0.05-1.png",sep="."), width=8, height=7, units="in", res=500)
density_plot <- ggplot(nam_snp_plot, aes(x=Divergence/1000000, color=Chr)) + 
  #geom_jitter(position=position_jitter(0.2),cex=0.005) +
  geom_density() + 
  #stat_peaks(colour = "red") +
  #stat_valleys(colour = "blue") + 
  #stat_summary(fun.y=median, geom="point", size=2) + 
  #geom_boxplot(width=0.1) +
  xlab("Divergence Time (million years)") + 
  #xlab("Line") + 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values=c("#87a2ba", "#334863", "#eda64a","#90766b","#ead09e","#ec853f","#00a950","#1c8b82","#ab6661","#79538e","#a2a5b4")) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) +
  #scale_x_continuous(limits = c(0,0.05),minor_breaks = seq(0,0.05, 0.005)) #+
  scale_x_continuous(limits = c(0.05,1),minor_breaks = seq(0.05,1, 0.01)) #+
#facet_wrap(~ Chr,scales="free", ncol =1, strip.position="right") 
density_plot_build <- ggplot_build(density_plot)

density_plot + stat_peaks(
  data = density_plot_build[['data']][[1]], # take a look at this object
  aes(x = x, y = density),
  colour = "red",
  size =1.5
) 
#nam_snp_plot_pericentromere<-nam_snp_plot_pericentromere[nam_snp_plot_pericentromere$Divergence >50000,]
dev.off()


png(paste(prefix,"wholechr.density_distribution.0-0.05.png",sep="."), width=8, height=7, units="in", res=500)
density_plot <- ggplot(nam_snp_plot, aes(x=Divergence/1000000, color=Chr)) + 
  #geom_jitter(position=position_jitter(0.2),cex=0.005) +
  geom_density() + 
  #stat_peaks(colour = "red") +
  #stat_valleys(colour = "blue") + 
  #stat_summary(fun.y=median, geom="point", size=2) + 
  #geom_boxplot(width=0.1) +
  xlab("Divergence Time (million years)") + 
  #xlab("Line") + 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values=c("#87a2ba", "#334863", "#eda64a","#90766b","#ead09e","#ec853f","#00a950","#1c8b82","#ab6661","#79538e","#a2a5b4")) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) +
  scale_x_continuous(limits = c(0,0.05),minor_breaks = seq(0,0.05, 0.001)) #+
#scale_x_continuous(limits = c(0.05,1),minor_breaks = seq(0.05,1, 0.05)) #+
#facet_wrap(~ Chr,scales="free", ncol =1, strip.position="right") 
density_plot_build <- ggplot_build(density_plot)

density_plot + stat_peaks(
  data = density_plot_build[['data']][[1]], # take a look at this object
  aes(x = x, y = density),
  colour = "red",
  size =1.5
) 
#nam_snp_plot_pericentromere<-nam_snp_plot_pericentromere[nam_snp_plot_pericentromere$Divergence >50000,]
dev.off()


png(paste(prefix,"wholechr.density_distribution.0-1.png",sep="."), width=10, height=7, units="in", res=500)
density_plot <- ggplot(nam_snp_plot, aes(x=Divergence/1000000, color=Chr)) + 
  #geom_jitter(position=position_jitter(0.2),cex=0.005) +
  geom_density() + 
  #stat_peaks(colour = "red") +
  #stat_valleys(colour = "blue") + 
  #stat_summary(fun.y=median, geom="point", size=2) + 
  #geom_boxplot(width=0.1) +
  xlab("Divergence Time (million years)") + 
  #xlab("Line") + 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values=c("#87a2ba", "#334863", "#eda64a","#90766b","#ead09e","#ec853f","#00a950","#1c8b82","#ab6661","#79538e","#a2a5b4")) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) +
  scale_x_continuous(limits = c(0,1),minor_breaks = seq(0,1, 0.01)) #+
  #facet_zoom(xy = Divergence > 0.05 & Divergence < 0.5) 
#scale_x_continuous(limits = c(0.05,1),minor_breaks = seq(0.05,1, 0.05)) #+
#facet_wrap(~ Chr,scales="free", ncol =1, strip.position="right") 
density_plot_build <- ggplot_build(density_plot)

density_plot + stat_peaks(
  data = density_plot_build[['data']][[1]], # take a look at this object
  aes(x = x, y = density),
  colour = "red",
  size =1.5
) 
#nam_snp_plot_pericentromere<-nam_snp_plot_pericentromere[nam_snp_plot_pericentromere$Divergence >50000,]
dev.off()




