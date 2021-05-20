library(karyoploteR)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(rstatix)

pangene_file<-'/scratch/jl03308/NAM_pancentromere/gene_synteny/pangene/pan-gene-expression-counts-tpm-rpkm.csv'
pangene<-read.csv(pangene_file)
NAMlines<-grep('count', colnames(pangene), value=TRUE) %>% str_replace("_count", "")
pangene_reformat <- data.frame()
for(line in NAMlines){
  x<-which(NAMlines==line)
  start_col<-(x-1)*9 + 2
  end_col<-x*9 +1
  tmp<-pangene[c(1,start_col:end_col)]
  colnames(tmp) <- c('PanGeneID','geneID','Chr','Start','End','Strand','Length','count', 'RPKM', 'TPM')
  tmp$line <- line
  pangene_reformat <- rbind(pangene_reformat,tmp)
}

pangene_reformat <- pangene_reformat %>% drop_na()

tpm_file<-'/scratch/jl03308/NAM_Canu1.8_verified_version_1/NAM_counts_tpm.txt'
tpm<-read.csv(tpm_file,sep ='\t', header = FALSE,col.names=c('line','geneID','Chr','Start','End','Strand','Length','TPM_new'))
#colnames(tmp) <- c('line','geneID','Chr','Start','End','Strand','Length','TPM_new')
pangene_merged<-merge(pangene_reformat, tpm, by=c("line",'geneID','Chr','Start','End','Strand','Length'))

umr_file<-'/scratch/jl03308/NAM_pancentromere/methylation/UMRs/UMRs_classification/NAM.line_UMR.TSSrange.bed'
umr<-read.csv(umr_file,sep ='\t', header = FALSE,col.names=c('line','Chr','UMR_Start','UMR_End','Methyl','Div','geneID','Strand'))
#colnames(tmp) <- c('line','geneID','Chr','Start','End','Strand','Length','TPM_new')
pangene_umr_merged<-merge(pangene_merged, umr, by=c("line",'geneID','Chr'))
pangene_merged_ref<-pangene_merged[pangene_merged$line=="B73",]
pangene_umr_merged_ref<-merge(pangene_umr_merged, pangene_merged_ref, by=c('PanGeneID','Chr'))
pangene_umr_merged_ref$log2FC <- abs(log2(pangene_umr_merged_ref$TPM_new.x + 0.01) -log2(pangene_umr_merged_ref$TPM_new.y+0.01))
pangene_umr_merged_ref_plot <- pangene_umr_merged_ref[,c('PanGeneID','Chr',"line.x",'Div','log2FC')]

col_vector=c("#D14799", "#A5AEDA", "#A17F89","#B88FD9","#9139EB","#E2859F","#DBBD4D","#BB66DF","#E5BBE1","#7FE5E4","#779A83","#C9D4E3","#7CE6BC","#69BCDE","#616BE0","#E3C3AD","#6590D6","#E845D4","#77E54E","#E25652","#D3EB58","#D2E9D5","#98CD73","#D0915A","#67EC94","#EB8CDE")

#pangene_umr_merged_ref[, c('Div')] <- sapply(pangene_umr_merged_ref[, c('Div')], as.numeric)
#pangene_umr_merged_ref[pangene_umr_merged_ref$TPM_diff >=200, ]$TPM_diff <- 200
pangene_umr_merged_ref_grouped <- pangene_umr_merged_ref_plot
pangene_umr_merged_ref_grouped$group[as.numeric(pangene_umr_merged_ref_grouped$Div) <=200000] <- '<=200K'
pangene_umr_merged_ref_grouped$group[as.numeric(pangene_umr_merged_ref_grouped$Div) >200000] <- '>200K'
pangene_umr_merged_ref_grouped$group[(pangene_umr_merged_ref_grouped$Div =="unaligned")] <- 'NeoUMR'

div<-c('<=200K','>200K','NeoUMR')
pangene_umr_merged_ref_grouped$group <- factor(pangene_umr_merged_ref_grouped$group,
                       levels = div)
pangene_umr_merged_ref_grouped<-unique(pangene_umr_merged_ref_grouped)
# for(chrs in chr){
#   pangene_umr_merged_ref_grouped_selected<-pangene_umr_merged_ref_grouped[(pangene_umr_merged_ref_grouped$Chr ==chrs) & ((pangene_umr_merged_ref_grouped$group =='NeoUMR') | (pangene_umr_merged_ref_grouped$group =='0-50K')),]
#   print(wilcox.test(log2FC ~ group, data = pangene_umr_merged_ref_grouped_selected, paired = FALSE))
# }
# 
# chr<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10')
# pangene_umr_merged_ref_grouped.test$Chr <- factor(pangene_umr_merged_ref_grouped.test$Chr,
#                        levels = chr)

library(plyr)
meds <- ddply(pangene_umr_merged_ref_grouped, .(group), summarize, med = median(log2FC))

options(bitmapType='cairo')

png('/scratch/jl03308/NAM_pancentromere/methylation/UMRs/UMRs_classification/test.png', width=8, height=6, units="in", res=500)
# ggplot(pangene_umr_merged_ref_grouped, aes(x=Div/1000000, y=TPM_diff, color=line.x)) +
#   geom_point(size = 0.00001) + 
ggplot(pangene_umr_merged_ref_grouped, aes(x=group, y=log2FC,color=group)) + 
  #geom_jitter(size=0.01,alpha=0.5) +
  #geom_violin(trim=FALSE) + 
  geom_boxplot(alpha=0) + 
  #stat_summary(fun.y = mean, geom = "errorbar", width = .75,color='red') + 
  scale_color_manual(values=c("#87a2ba", "#AB6661", "#eda64a","#666666","#ead09e","#604444","#2f4335","#1c8b82","#ab6661","#79538e","#a2a5b4","#dfc064","#90766b","#f4a3a4","#ffe2e2","#ec853f","#60965b")) +
  xlab("UMR types (NAM/B73)") + 
  ylab("log2(FC) of NAM/B73") + 
  #scale_y_continuous(trans='log10') + 
  #scale_color_manual(values=col_vector) +
  #theme(legend.position="top") +
  geom_text(data = meds, aes(y = med, label = round(med,2)),size = 5) + 
  theme_minimal() +
  theme(axis.text.x=element_text(size=13),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=16,face="bold"), 
        strip.text.x = element_text(size =15, colour = "black"),
        strip.text.y = element_text(size =15, colour = "black"),
        legend.position = "none") +
  #scale_x_continuous(limits = c(0,1)) + 
  scale_y_continuous(limits = c(0,2.5)) #+ 
  #facet_grid(Chr~.) 
dev.off()

