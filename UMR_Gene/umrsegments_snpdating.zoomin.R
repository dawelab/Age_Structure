#!/usr/bin/env Rscript
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
wd_dir <- args[1]
gff_file <- args[2]
centromere_file <- args[3]
genome_size_file <- args[4]
umr_file <- args[5]
chr <- args[6]
ref <- args[7]
group <- args[8]
start_coord<- args[9]
end_coord<- args[10]
div_file<- args[11]

start_coord<- as.numeric(start_coord)
end_coord<- as.numeric(end_coord)

library(karyoploteR)
library(data.table)
library(ggplot2)
library(dplyr) 
library(stringr)
library(tidyr)
options(bitmapType='cairo')

setwd(wd_dir)
pp <- getDefaultPlotParams(plot.type = 1)
#data1height 400 for NAM group
pp$data1height <- 400
pp$ideogramheight <- 0.5
pp$data1inmargin<-0
topmargin=0
bottommargin=0
data1outmargin=0
# wd_dir <- '/scratch/jl03308/NAM_pancentromere/NAM_SV/chr9'
# gff_file <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/annotation/pan_TE_v1.0/B73.PLATINUM.pseudomolecules-v1.fasta.mod.EDTA.TEanno.gff'
# chr<-'chr9'
# centromere_file <- '/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum'
# genome_size_file <- '/scratch/jl03308/NAM_Canu1.8_verified_version_1/pseudomolecules/B73.PLATINUM.pseudomolecules-v1.chrom.sizes'
# ref<-'B73'
# alignment_file <- '/scratch/jl03308/NAM_pancentromere/NAM_SV/chr2/B73_NAM.aligned.bed'

gff<-read.table(gff_file)
knob180 <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="knob/knob180"),c('V1','V4','V5')]
TR1 <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="knob/TR-1"),c('V1','V4','V5')]
subtelomere <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="subtelomere/4-12-1"),c('V1','V4','V5')]
CentC <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="Cent/CentC"),c('V1','V4','V5')]
nor <- gff[(sapply(strsplit(as.character(gff$V1), "_"), "[", 2)==chr) & (gff$V3=="rDNA/spacer"),c('V1','V4','V5')]

#centromere_file<-'/scratch/jl03308/NAM_pancentromere/analysis/peak_call/NAM/NAM.centro.coords.assembly.status.content.sum'
centromere <- read.table(centromere_file)
ref.centromere <- centromere[(centromere$V1==ref)&(centromere$V2==chr),c('V1','V3','V4')]

#TR1 <- read.table(TR1_file)
#'/Users/jianingliu/Downloads/pan_TE_v1.0/knob180.bed'
#knob180 <- read.table(knob180_file)
#subtelomere <- read.table(subtelomere_file)
#CentC <- read.table(centc_file)
#B73.chr.centromere <- data.frame("chr" = c('B73'), "Start" = c(50495000), "End" = c(52384999))

#TR1.chr<-TR1[TR1$V1==chr,]
if(nrow(TR1) > 0){TR1$V1 <- ref}
#knob180.chr<-knob180[knob180$V1==chr,]
if(nrow(knob180) > 0){knob180$V1 <-ref}
#subtelomere.chr<-subtelomere[subtelomere$V1==chr,]
if(nrow(subtelomere) > 0){subtelomere$V1 <-ref}
#CentC.chr<-CentC[CentC$V1==chr,]
if(nrow(CentC) > 0){CentC$V1 <-ref}
if(nrow(nor) > 0){nor$V1 <-ref}

div_file<-'/scratch/jl03308/NAM_pancentromere/NAM_SV/chr2/B73_NAM.div.txt'
div<-read.table(div_file)

#umr_file <- '/scratch/jl03308/NAM_pancentromere/methylation/UMRs/NAM.ref_B73.UMR.syn.bed'
umr.bed<-read.table(umr_file)
umr.bed<-umr.bed[(umr.bed$V1 ==chr),]
all_umr=umr.bed[,c('V1','V2','V3')]
all_umr$V1<-ref
all_umr=all_umr[(all_umr$V2>=start_coord)&(all_umr$V3<=end_coord),]

# info = file.info(core_file)
# if(info$size ==0){core<-data.frame()}else{
#   core <- read.table(core_file)
#   core$V1<-'B73'
# }

pangene_file<-'/scratch/jl03308/NAM_pancentromere/gene_synteny/pangene/pan-gene-expression-counts-tpm-rpkm.csv'
pangene<-read.csv(pangene_file)
#head(pangene,2)
#pangene<-pangene[-c(1)]

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
  #print(line)
  #print(x)
}
# identify all B73 genes (present), and the corresponding pangene ID
refgenes <-pangene_reformat[(pangene_reformat$line=='B73')&(pangene_reformat$Chr==chr)&!is.na(pangene_reformat$Start),c('line','Start','End','Strand')]
refgenes <- refgenes[(refgenes$Start>=start_coord)&(refgenes$End<=end_coord),]
refgenes_forward <- refgenes[(refgenes$Strand=="+"),c('line','Start','End')]
refgenes_reverse <- refgenes[(refgenes$Strand=="-"),c('line','Start','End')]

#genome_size<-read.table(genome_size_file)
genome<-data.frame("line" = c('B73'), "Start" = c(start_coord), "End" = c(end_coord))

png(paste(chr,group,start_coord, end_coord, "SNPdating.UMR.png",sep="."), width=10, height=8, units="in", res=500)
#pdf(paste(chr,group,start_coord, end_coord, "SNPdating.UMR.pdf",sep="."), width = 10, height = 6)

chromosome <- plotKaryotype(genome = genome, plot.type=1, plot.params = pp)
#start=0.025
#end=0.05

#kpRect(chromosome, data=toGRanges(all_umr),x0=all_umr$V2,x1=all_umr$V3,y0=0, y1=1,lwd=0.00,r0=start, r1=end,col='black',border=NA,clipping = TRUE)
#kpAddLabels(chromosome, labels='Total', r0=start, r1=end,side="left",label.margin=0.005,cex=.9)

kpRect(chromosome, data=toGRanges(refgenes_forward),x0=refgenes_forward$Start,x1=refgenes_forward$End,y0=0, y1=1,lwd=0.00,r0=0, r1=0.04,col='#aa165b',border=NA)
kpRect(chromosome, data=toGRanges(refgenes_reverse),x0=refgenes_reverse$Start,x1=refgenes_reverse$End,y0=0, y1=1,lwd=0.00,r0=0, r1=0.04,col='#0f7646',border=NA)
kpAddBaseNumbers(chromosome, tick.dist = 1000000, add.units = FALSE,digits=1,cex=0,tick.len = 8,minor.tick.dist = 100000, minor.tick.len = 3,clipping=TRUE)

start=0.04
end=0.065
query=unique(umr.bed$V4)

for(x in query){
  line_specific_umr=umr.bed[(umr.bed$V4==x)&(umr.bed$V2>=start_coord)&(umr.bed$V3<=end_coord),c('V4','V2','V3')]
  line_specific_umr$V4<-ref
  kpRect(chromosome, data=toGRanges(line_specific_umr),x0=line_specific_umr$V2,x1=line_specific_umr$V3,y0=0, y1=1,lwd=0.00,r0=start, r1=end,col='black',border=NA,clipping = TRUE)
  kpAddLabels(chromosome, labels=x, r0=start-0.01, r1=end-0.01,side="left",label.margin=0.005,cex=.9)
  start=start+0.05
  end=end+0.05
}

start=0.065
end=0.09
query=unique(div$V2)
for(x in query){
  #print(x)
  #x<-'B97'
  line_specific=div[(div$V2==x)&(div$V3>=start_coord)&(div$V4<=end_coord),c('V1','V3','V4','V7')]
  line_specific[line_specific$V7>600000,]$V7=600000
  #quantile(div[!is.na(div$V5),]$V5, probs = 0.02,names = FALSE)
  kpHeatmap(chromosome, toGRanges(line_specific), y=line_specific$V7, color=c("#e7eaed",'#e6ab00',"#fb7c24",'#cc0000','#6d1c49',"#522b6c","black"), ymin=0, ymax=600000,r0=start,r1=end,clipping = TRUE)
  #kpAddLabels(chromosome, labels=x, r0=start, r1=end,side="left",label.margin=0.005,cex=.9)
  start=start+0.05
  end=end+0.05
}
dev.off()
