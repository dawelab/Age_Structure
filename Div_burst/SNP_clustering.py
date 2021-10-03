#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from itertools import permutations, islice
import random
import pandas as pd
import numpy as np
import seaborn as sns
import prince
import re
import os
from sklearn.preprocessing import normalize
import scipy.cluster.hierarchy as shc
import matplotlib
import matplotlib.pyplot as plt

colors=["#334863", "#eda64a","#90766b","#ead09e","#ec853f","#00a950","#1c8b82","#ab6661","#79538e","#a2a5b4"]
chrs=['chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10', 'chr1']
matplotlib.rcParams['lines.linewidth'] = 3.5

for i in range(len(chrs)):
    chr=chrs[i]
    color=colors[i]
    print(chr)
    synsnp_file='/scratch/jl03308/NAM_pancentromere/pairwise_NAM_SV/' + chr + '/B73.NAM.pericentromere.syn.intergenic.allsnps.sum'
    ref=os.path.basename(synsnp_file).split('.')[0]
    syn_snps=pd.read_csv(synsnp_file, delimiter="\t",comment="#",header=None).sort_values([2,0],ascending=[True,True]).drop_duplicates().reset_index(drop=True)
    syn_snps_reshape=pd.DataFrame()
    #row index is NAM line, column index is SNP coords
    snp_coord_reshape=syn_snps.pivot(index=0, columns=2,values=5)
    snp_coord_reshape_pca=snp_coord_reshape.replace("T", 1).replace("A", 2).replace("C", 3).replace("G", 4)
    snp_coord_reshape_pca=snp_coord_reshape_pca.dropna(axis='columns')
    snp_coord_reshape_scaled = normalize(snp_coord_reshape_pca)
    snp_coord_reshape_scaled = pd.DataFrame(snp_coord_reshape_scaled, columns=snp_coord_reshape_pca.columns,index=snp_coord_reshape_pca.index)
    outdir=os.path.dirname(synsnp_file)
    prefix='.'.join(os.path.basename(synsnp_file).split('.')[0:5])
    figfile=os.path.join(outdir, '.'.join((prefix, chr, "dendrogram.png")))
    #snp_coord_reshape_allchrs=PCA_analysis(snp_coord_reshape_allchrs,figfile)
    plt.figure(figsize=(8, 10))
    #dend = shc.dendrogram(shc.linkage(snp_coord_reshape_scaled, method='ward'),labels=snp_coord_reshape_scaled.index,color_threshold=0, above_threshold_color="#383637")
    dend = shc.dendrogram(shc.linkage(snp_coord_reshape_scaled, method='ward'),labels=snp_coord_reshape_scaled.index,color_threshold=0, above_threshold_color=color, leaf_font_size=13)
    plt.savefig(figfile,markerfacecolor="None",transparent=True)
