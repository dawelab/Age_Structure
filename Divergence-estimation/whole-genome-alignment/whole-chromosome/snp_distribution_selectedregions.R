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
pericentromere<-rbind(c('chr1',125000000,135000000),
                      c('chr4',110000000,120000000),
                      c('chr5',102000000,112000000),
                      c('chr9',58000000,65000000),
                      c('chr10',47000000,52000000))

pericentromere<-data.frame(pericentromere)
colnames(pericentromere) <- c('Chr','Start','End')

nam_snp_plot_pericentromere<-data.frame()
for(x in 1:nrow(pericentromere)){
  #x<-1
  chr_selected <- pericentromere[x,'Chr']
  chr_selected_start <- as.numeric(pericentromere[x,'Start'])
  chr_selected_end <- as.numeric(pericentromere[x,'End'])
  nam_snp_plot_chr_selected<-nam_snp_plot[(nam_snp_plot$Chr==chr_selected) & (nam_snp_plot$Start>=chr_selected_start) & (nam_snp_plot$End<=chr_selected_end),]
  nam_snp_plot_pericentromere<-rbind(nam_snp_plot_pericentromere,nam_snp_plot_chr_selected)
}

#output_file<-'NAM.intergeni_window_div'
# png(paste(output_file,"pericentromere.distribution.png",sep="."), width=8, height=15, units="in", res=500)
# #nam_snp_plot_pericentromere<-nam_snp_plot_pericentromere[nam_snp_plot_pericentromere$Divergence >50000,]
# ggplot(nam_snp_plot_pericentromere, aes(x=NAM_line, y=Divergence/1000000)) + 
#   #geom_jitter(position=position_jitter(0.2),cex=0.005) +
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

#nam_snp <- data.frame(read.table(snp_file))
png(paste(prefix,".selectedregions.pericentromere.png",sep="."), width=15, height=22, units="in", res=500)
ggplot(nam_snp_plot_pericentromere, aes(x=Start/1000000, y=Divergence/1000000, color=NAM_line)) +
  geom_point(cex = 0.3,alpha=0.5) + 
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
  scale_y_continuous(limits = c(0,0.5)) + 
  facet_wrap(~ Chr,scales="free", ncol =1, strip.position="right") 
#facet_grid(Chr~.,scales="free")
dev.off()

