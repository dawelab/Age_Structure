#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)
library(viridis)
#install.packages("ggpubr")
library(ggpubr)
options(bitmapType='cairo')

wdir <-'/scratch/jl03308/NAM_pancentromere/methylation/UMRs/UMRs_classification'
setwd(wdir)
umr_file<-'NAM.line_UMR.classification.bed'

umr <- unique(data.frame(read.table(umr_file)))
#umr$V2 <- sub("^", "chr", umr$V2 )
colnames(umr) <- c("NAM_line",'Chr','Start','End','Methyl','Div')
umr$Type[umr$Div >200000] <-'Highdiv (>200K)'
umr$Type[umr$Div <=200000] <-'Lowdiv (<=200K)'
umr$Type[umr$Div=="unaligned"] <-'Neo'
umr$Size <-abs(umr$End-umr$Start) +1

umr_plot <- umr[,c("NAM_line",'Chr','Lowdiv','Highdiv (>200K)','Neo')]
chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
umr_plot<-melt(umr_plot, id=c("NAM_line", "Chr"))
umr_plot$Chr <- factor(umr_plot$Chr,
                       levels = chr)

#nam_centromere_plot$Chr <- factor(nam_centromere_plot$Chr,levels = chr)
#png('NAM.UMR_classification.png', width = 8, height = 8, units = "in",res=300)
#pdf hide white lines in between stacked bar plots
pdf('NAM.UMR_classification.pdf', width = 8, height = 8)
ggplot(umr, aes(fill=Type,x=NAM_line, y=Size/1000000)) + 
  geom_bar(position="stack", stat="identity",colour =NA, size = 0) + 
  #stat_summary(fun.y=median, geom="point", size=2) + 
  ylab("UMR Size (Mb)") + 
  xlab("Line") + 
  #scale_y_continuous(trans='log10') + 
  scale_fill_manual(values=c("#87a2ba", "#AB6661", "#eda64a","#666666","#ead09e","#604444","#2f4335","#1c8b82","#ab6661","#79538e","#a2a5b4","#dfc064","#90766b","#f4a3a4","#ffe2e2","#ec853f","#60965b")) +
  labs(fill = "") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) #+
  #facet_grid(Chr~.) 
dev.off()

wdir <-'/scratch/jl03308/NAM_pancentromere/methylation/UMRs/UMRs_classification'
setwd(wdir)
umr_file<-paste('NAM.UMR.classified.bed',sep="")

umr <- data.frame(read.table(umr_file))
aligned_umr<- umr[umr$V6 != "unaligned",]
unaligned_umr <- umr[umr$V6 == "unaligned",]
colnames(aligned_umr) <- c("NAM_line",'Chr','Start','End','Methyl','Div')
#aligned_umr_plot <- umr[,c("NAM_line",'Chr','Div')]
chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
aligned_umr$Chr <- factor(aligned_umr$Chr,
                       levels = chr)
aligned_umr[, c('Div')] <- sapply(aligned_umr[, c('Div')], as.numeric)
#aligned_umr$div<- as.numeric(as.character(aligned_umr$div))
pdf('NAM.alignedUMR.div2.pdf', width = 8, height = 12)

library(RColorBrewer)
n <- 27
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#png('NAM.alignedUMR.div.png', width = 8, height = 12, units = "in",res=300)
# ggplot(aligned_umr, aes(x=NAM_line, y=Div/1000000)) + 
#   ggplot(aligned_umr, aes(x=Div/1000000, color=NAM_line)) + 
#   geom_density() + 
#   geom_histogram(aes(y=..density..), alpha=0.5, 
#                  position="identity") + 
#   #stat_summary(fun.y=median, geom="point", size=2) + 
#   ylab("Divergence (million years)") + 
#   xlab("Line") + 
#   #scale_y_continuous(trans='log10') + 
#   scale_fill_manual(values=col_vector) +
#   labs(fill = "") +
#   theme_minimal() +
#   theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
#         axis.text.y=element_text(size=13),
#         axis.title=element_text(size=16,face="bold"), 
#         strip.text.x = element_text(size =15, colour = "black"),
#         strip.text.y = element_text(size =15, colour = "black")) +
#   facet_grid(Chr~.) 
# dev.off()

col_vector=c("#D14799", "#A5AEDA", "#A17F89","#B88FD9","#9139EB","#E2859F","#DBBD4D","#BB66DF","#E5BBE1","#7FE5E4","#779A83","#C9D4E3","#7CE6BC","#69BCDE","#616BE0","#E3C3AD","#6590D6","#E845D4","#77E54E","#E25652","#D3EB58","#D2E9D5","#98CD73","#D0915A","#67EC94","#EB8CDE")
png('NAM.alignedUMR.div.png', width = 8, height = 12, units = "in",res=300)
ggplot(aligned_umr, aes(x=Div/1000000, color=NAM_line)) + 
  stat_ecdf(geom = "step",size=0.3) + 
  #stat_summary(fun.y=median, geom="point", size=2) + 
  ylab("Percentage of UMR") + 
  xlab("Divergence (million years)") + 
  #scale_y_continuous(trans='log10') + 
  scale_color_manual(values=col_vector) +
  labs(fill = "") +
  theme_minimal() +
  theme(axis.text.x=element_text(size=13,angle = 40, hjust = 1),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black")) +
  facet_grid(Chr~.) 
dev.off()


# library(randomcoloR)
# n <- 27
# palette <- distinctColorPalette(n)
