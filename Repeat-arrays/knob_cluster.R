#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
wdir <- args[1]
pairwise_distance_file <- args[2]
output_file <- args[3]

library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(reshape2)
library(tools)
options(bitmapType='cairo')

#wdir<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/'
setwd(wdir)
#pairwise_distance_file<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/NAM.knob180+TR-1.chr9.knob1.pairwise_distance.txt'
#output_file<-'NAM.knob180.chr3.knob1'
pairwise_distance <- data.frame(read.table(pairwise_distance_file))
pairwise_distance <- transform(pairwise_distance, minlen = pmin(V6, V7))
pairwise_distance$V11 <- pairwise_distance$V8/pairwise_distance$minlen
# lines <- sort(unique(unlist(pairwise_distance[,c('V4','V5')])))
# pairwise_distance_square <- matrix(NA, length(lines), length(lines), dimnames = list(lines, lines))
# pairwise_distance_square[as.matrix(pairwise_distance[c('V4','V5')])] <- pairwise_distance[["V11"]]
# pairwise_distance_square[as.matrix(pairwise_distance[c('V5','V4')])] <- pairwise_distance[["V11"]]
# pairwise_distance_square[is.na(pairwise_distance_square)] <- 0

#pairwise_distance_matrix <- acast(pairwise_distance[,c('V4','V5','V11')], V5 ~ V4,value.var = "V11")
#pairwise_distance_matrix[is.na(pairwise_distance_matrix)] <- 0
#pairwise_dist <- as.dist(pairwise_distance_square, diag = TRUE)

type<-unlist(strsplit(file_path_sans_ext(basename(pairwise_distance_file)),"\\."))[[2]]
prefix<-unlist(strsplit(file_path_sans_ext(basename(pairwise_distance_file)),"\\."))[[4]]
chr<-unlist(strsplit(file_path_sans_ext(basename(pairwise_distance_file)),"\\."))[[3]]

repeat_num_file<-'/scratch/jl03308/NAM_pancentromere/analysis/syn_repeat/dating/NAM.repeat.num.txt'
repeat_num <- data.frame(read.table(repeat_num_file))
if(type=="knob180+TR-1"){
  repeat_num_selected <- repeat_num[((repeat_num$V1=="knob180") |(repeat_num$V1=="TR-1")) &(repeat_num$V2==chr)&(repeat_num$V3==prefix),]
  } else{repeat_num_selected <- repeat_num[(repeat_num$V1==type)&(repeat_num$V2==chr)&(repeat_num$V3==prefix),]
}
nam_lines <- unique(repeat_num_selected$V4)
pairwise_distance_square <- matrix(NA, length(nam_lines), length(nam_lines), dimnames = list(nam_lines, nam_lines))
pairwise_distance_square[as.matrix(pairwise_distance[c('V4','V5')])] <- pairwise_distance[["V11"]]
pairwise_distance_square[as.matrix(pairwise_distance[c('V5','V4')])] <- pairwise_distance[["V11"]]
pairwise_distance_square[is.na(pairwise_distance_square)] <- 0

temperate_origin=c('B73','B97','Ky21','M162W','MS71','Oh43','Oh7b')
mixed_origin = c('M37W','Mo18W','Tx303')
flint_origin=c('HP301','P39','IL14H')
tropical_origin = c('CML52','CML69','CML103','CML228','CML247','CML277','CML322','CML333','Ki3','Ki11','NC350','NC358','Tzi8')

groups<-list(absence = which(nam_lines %in% repeat_num_selected[repeat_num_selected$V5==0,]$V4),
     temperate=which(nam_lines %in% temperate_origin[temperate_origin %in% repeat_num_selected[repeat_num_selected$V5>0,]$V4]),
     mixed=which(nam_lines %in% mixed_origin[mixed_origin %in% repeat_num_selected[repeat_num_selected$V5>0,]$V4]),
     flint=which(nam_lines %in% flint_origin[flint_origin %in% repeat_num_selected[repeat_num_selected$V5>0,]$V4]),
     tropical=which(nam_lines %in% tropical_origin[tropical_origin %in% repeat_num_selected[repeat_num_selected$V5>0,]$V4]))

#colnames(nam_centromere) <- c("NAM_line", "Chr","Start","End","Size","100N gap number","13N gap number","Known gap size","Fully assembled (Y/n)","CentC","CRM2","CRM1","cinful-zeon","opie-ji","huck","prem1","grande", "All_TE","gene","Known sequences", "Other TE","Total gap","Unknown" )
# hc1 <- hclust(pairwise_dist)
# png(paste(output_file,"hc.png",sep="."), width=5, height=4, units="in", res=500)
# # Plot the obtained dendrogram
# plot(hc1, cex = 0.6, hang = -1)
# #nam_centromere_plot$Chr <- factor(nam_centromere_plot$Chr,levels = chr)
# dev.off()
#got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

# groups <- list(temperate = c('B73','B97','Ky21','M162W','MS71','Oh43','Oh7b'),
#                mixed = c('M37W','Mo18W','Tx303'),
#                flint=c('HP301','P39','IL14H'),
#                tropical = c('CML52','CML69','CML103','CML228','CML247','CML277','CML322','CML333','Ki3','Ki11','NC350','NC358','Tzi8'))
# groups <- list(temperate = c(1,2,15,16,19,22,23),
#                mixed = c(17,18,25),
#                flint=c(11,12,24),
#                tropical = c(3,4,5,6,7,8,9,10,13,14,20,21,26))

library(qgraph)

png(paste(output_file,"network.png",sep="."), width=5, height=5, units="in", res=500)
qgraph(pairwise_distance_square, layout='spring', vsize=4, labels = colnames(pairwise_distance_square),
       groups = groups,
       color = c("#eee9ee", "#58B5BC", "#F5651C","#E99093",'#93a7da'),
       edge.color = "darkgray", edge.width = 0.8,legend=FALSE) 
    
# ggraph(pairwise_dist,layout = "spring") + 
#   theme_graph()+
#   theme(legend.position = "none")  
  
dev.off()
