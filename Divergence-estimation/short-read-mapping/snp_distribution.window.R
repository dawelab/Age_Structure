#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
wdir <- args[1]
snp_file <- args[2]
output_file <- args[3]

library(ggplot2)
options(bitmapType='cairo')

setwd(wdir)

nam_snp <- data.frame(read.table(snp_file))
nam_snp_plot <- nam_snp[,c('V1','V2','V5','V7')]
#colnames(nam_centromere) <- c("NAM_line", "Chr","Start","End","Size","100N gap number","13N gap number","Known gap size","Fully assembled (Y/n)","CentC","CRM2","CRM1","cinful-zeon","opie-ji","huck","prem1","grande", "All_TE","gene","Known sequences", "Other TE","Total gap","Unknown" )
colnames(nam_snp_plot) <- c("NAM_line", 'Chr',"SNP",'Divergence')
chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
nam_snp_plot$Chr <- factor(nam_snp_plot$Chr,
                           levels = chr)

#nam_centromere_plot$Chr <- factor(nam_centromere_plot$Chr,levels = chr)
png(paste(output_file,"distribution.png",sep="."), width=8, height=15, units="in", res=500)

ggplot(nam_snp_plot, aes(x=NAM_line, y=Divergence/1000000)) + 
  geom_violin() + 
  #stat_summary(fun.y=median, geom="point", size=2) + 
  #geom_boxplot(width=0.1) +
  ylab("Divergence Time (million years)") + 
  xlab("Line") + 
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
  facet_grid(Chr~.) 
dev.off()

library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

#nam_snp <- data.frame(read.table(snp_file))
nam_snp_plot <- nam_snp[,c('V1','V2','V4','V7')]
colnames(nam_snp_plot) <- c("NAM_line", 'Chr',"Loc",'Divergence')
chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
nam_snp_plot$Chr <- factor(nam_snp_plot$Chr,
                           levels = chr)
png(paste(output_file,"wholechr.png",sep="."), width=15, height=15, units="in", res=500)
ggplot(nam_snp_plot, aes(x=Loc/1000000, y=Divergence/1000000, color=NAM_line)) +
  geom_point(size = 0.00001) + 
  ylab("Divergence Time (million years)") + 
  xlab("Chromosome (Mb)") + 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values=col_vector) +
  #theme(legend.position="top") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) +
  scale_y_continuous(limits = c(0,1)) + 
  facet_grid(Chr~.) 
dev.off()

# nam_snp_chr<-nam_snp_plot[(nam_snp_plot$Chr=="chr6")&(nam_snp_plot$Loc>=15000000)&(nam_snp_plot$Loc<=25000000),]
# png(output_file, width=15, height=5, units="in", res=500)
# ggplot(nam_snp_chr, aes(x=Loc/1000000, y=Divergence/1000000, color=NAM_line)) +
#   geom_point(size = 0.1) + geom_line() + 
#   ylab("Divergence Time (million years)") + 
#   xlab("Mb") + 
#   #scale_y_continuous(trans='log10') + 
#   scale_color_manual(values=col_vector) +
#   #theme(legend.position="top") +
#   theme_minimal() +
#   theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
#         axis.text.y=element_text(size=13),
#         axis.title=element_text(size=16,face="bold"), 
#         strip.text.x = element_text(size =15, colour = "black"),
#         strip.text.y = element_text(size =15, colour = "black")) +
#   scale_y_continuous(limits = c(0,1)) + 
#   facet_grid(NAM_line~.) 
# dev.off()
