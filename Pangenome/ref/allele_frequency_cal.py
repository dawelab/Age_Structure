#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#from itertools import permutations, islice
import random
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
#import re
import os


#NAMline=['B73','B97','CML103','CML228','CML247','CML277','CML322','CML333','CML52','CML69','HP301','IL14H','Ki11','Ki3','Ky21','M162W','M37W','Mo18W','MS71','NC350','NC358','Oh43','Oh7b','P39','Tx303','Tzi8']
flint_origin=['HP301','P39','IL14H']
temperate_origin=['B97','Ky21','M162W','MS71','Oh43','Oh7b']
mixed_origin = ['M37W','Mo18W','Tx303']
tropical_origin = ['CML52','CML69','CML103','CML228','CML247','CML277','CML322','CML333','Ki3','Ki11','NC350','NC358','Tzi8']
NAMline_origin=flint_origin+temperate_origin+mixed_origin+tropical_origin

def shuffled(array):
    random.shuffle(array)
    return list(array)

def main(alignment_file,genome_file,ref,outdir,frequency_file,saturation_stat_file):
    # read syn alignment
    #alignment_file='/Users/jianingliu/Downloads/B73.syn.aligned.bed'
    #ref='B73'
    #genome_file='/Users/jianingliu/Downloads/NAM.genome'
    genome=pd.read_csv(genome_file, delimiter="\t",comment="#",header=None)
    interval_pangenome_allchrs=pd.DataFrame()
    aligned_segments=pd.read_csv(alignment_file, delimiter="\t",comment="#",header=None).sort_values([10,1,2],ascending=[True,True,False]).reset_index(drop=True)
    chrs=aligned_segments[10].unique()
    #chrs=['chr10','chr9']
    NAMline=aligned_segments[3].unique()
    NAMline.sort()
    for selected_chr in chrs:
        #selected_chr='chr9'
        print("Detecting breakpoints " + selected_chr)
        aligned_segments_chr=aligned_segments[aligned_segments[10]==selected_chr]
        # extract start breakpoints
        breakpoint_start=aligned_segments_chr[[1,3]]
        # extract end breakpoints
        breakpoint_end=aligned_segments_chr[[2,3]]
        # merge start and end breakpoints
        breakpoint_end=breakpoint_end.rename({2: 1}, axis='columns')
        breakpoints=breakpoint_start.append(breakpoint_end).drop_duplicates().reset_index(drop=True)
        # start and end coord of each chr is the final coord breakpoints
        #start=1
        #end=182411202
        #
        genome_selected=genome[(genome[0]==selected_chr)&(genome[1]==ref)][[1,2,3]]
        start=1
        end=genome_selected[3].values[0]
        breakpoints=np.unique(np.concatenate((aligned_segments_chr[1].unique(), aligned_segments_chr[2].unique())))
        breakpoints=np.sort(np.unique(np.append(breakpoints,[start,end])))

        '''if breakpoint interval is smaller than 20bp, get the mid-point'''
        adjusted=[]
        adjusted_index=[]
        breakpoint_interval=np.diff(breakpoints)
        #breakpoints[:20]
        for i in range(0,len(breakpoint_interval)):
            if breakpoint_interval[i] < 20:
                adjusted_index += [i,i+1]
                adjusted += [int((breakpoints[i] + breakpoints[i+1])/2)]
        #adjusted[:20]
        #adjusted_index[:20]
        breakpoint_adjusted=np.sort(np.unique(np.concatenate((np.delete(breakpoints, adjusted_index), np.array(adjusted)))))
        if breakpoint_adjusted[0]!=start and abs(breakpoint_adjusted[0]-start) < 20:
            breakpoint_adjusted[0]=1
        else:
            breakpoint_adjusted=np.concatenate((breakpoint_adjusted, np.array([start])))
        if breakpoint_adjusted[-1]!=end and abs(breakpoint_adjusted[-1]-end) < 20:
            breakpoint_adjusted[-1]=end
        else:
            breakpoint_adjusted=np.concatenate((breakpoint_adjusted, np.array([end])))

        print("Calculate alignment frequency " + selected_chr)
        '''derive interval from breakpoints'''
        breakpoint_adjusted=np.sort(np.unique(breakpoint_adjusted))
        interval=[]
        for i in range(1,len(breakpoint_adjusted)):
            interval.append([breakpoint_adjusted[i-1],breakpoint_adjusted[i]])
        #interval[:20]
        #breakpoint_adjusted[breakpoint_adjusted==272741]
        interval_pangenome=pd.DataFrame(interval)
        #pd.DataFrame(interval_pangenome[1]-interval_pangenome[0]).boxplot()

        '''print the aligned coords in query if alignable'''
        for x in NAMline:
            #x='P39'
            aligned_segments_line=aligned_segments_chr[aligned_segments_chr[3]==x]
            for j in aligned_segments_line.index:
                #print(j)
                #aligned_segments_line.loc[j]
                start_coord=aligned_segments_line.loc[j,1]
                end_coord=aligned_segments_line.loc[j,2]
                query_start=aligned_segments_line.loc[j,4]
                query_end=aligned_segments_line.loc[j,5]
                interval_pangenome_start_end=interval_pangenome[(interval_pangenome[0]>=start_coord-20)&(interval_pangenome[1]<=end_coord+20)]
                if len(interval_pangenome_start_end) >0:
                    interval_pangenome.loc[interval_pangenome_start_end.index,x]=str(query_start) + ":" + str(query_end)
        '''convert to binary format'''
        interval_pangenome_binary=interval_pangenome.fillna(0)
        for i in range(len(NAMline)):
            interval_pangenome_binary.loc[interval_pangenome_binary[NAMline[i]]!=0,NAMline[i]] =1
        interval_pangenome_binary['chr'] = selected_chr
        interval_pangenome_allchrs=interval_pangenome_allchrs.append(interval_pangenome_binary)

    '''summary stat'''
    '''Write interval frequency, input for bar plotting, and gene/UMR analsyis'''
    print("Writing alignment frequency" )
    interval_pangenome_allchrs['sum'] = interval_pangenome_allchrs[NAMline].sum(axis=1)
    # sum_stat=[]
    # for i in range(0,len(NAMline)+1):
    #     #i=0
    #     alignedlen=(interval_pangenome_binary[interval_pangenome_binary['sum']==i][1] - interval_pangenome_binary[interval_pangenome_binary['sum']==i][0] + 1).sum()
    #     sum_stat.append([i,alignedlen])
    #pd.DataFrame(sum_stat).plot.bar(x=0, y=1)
    #frequency_file='/Users/jianingliu/Downloads/conservation.sum.txt'
    interval_pangenome_allchrs.loc[:,[0,1,'sum']]= interval_pangenome_allchrs.loc[:,[0,1,'sum']].astype('int64')
    interval_pangenome_allchrs=interval_pangenome_allchrs.rename({0: 'start',1:'end'}, axis='columns')
    pd.DataFrame(interval_pangenome_allchrs).to_csv(frequency_file,index=False,header=True,sep='\t')

    print("Plotting pairwise distance")
    #sns.set_theme(style="white")
    #a = 5  # number of rows
    #b = 2  # number of columns
    #fig, ax = plt.subplots(nrows=a, ncols=b,figsize=(20, 100),facecolor='w')
    for selected_chr in interval_pangenome_allchrs['chr'].unique():
        #c = int(re.sub('chr', '', selected_chr)) -1
        '''creat similarity matrix based on alignment 0/1'''
        similarity_matrix=pd.DataFrame(index=NAMline_origin,columns=NAMline_origin)
        interval_pangenome_binary=interval_pangenome_allchrs[interval_pangenome_allchrs['chr']==selected_chr]
        for i in range(len(NAMline_origin)-1):
            #i=0
            for j in range(i+1,len(NAMline_origin)):
                #j=1
                shared=interval_pangenome_binary[(interval_pangenome_binary[NAMline_origin[i]]==1)&(interval_pangenome_binary[NAMline_origin[j]]==1)]
                sharedlen=(shared['end']-shared['start']+1).sum()
                total=interval_pangenome_binary[(interval_pangenome_binary[NAMline_origin[i]]==1)|(interval_pangenome_binary[NAMline_origin[j]]==1)]
                totallen=(total['end']-total['start']+1).sum()
                similarity=sharedlen/totallen
                similarity_matrix.loc[NAMline_origin[j],NAMline_origin[i]] = similarity

        similarity_matrix=similarity_matrix.fillna(0)
        similarity_matrix_plot=np.array(similarity_matrix)
        '''plot similarity matrix'''
        #similarity_matrix = similarity_matrix[NAMline_origin]
        #plt.subplot(a, b, c)

        mask = np.triu(np.ones_like(similarity_matrix_plot, dtype=bool))
        # Generate a custom diverging colormap
        cmap = sns.diverging_palette(230, 20, as_cmap=True)
        # Draw the heatmap with the mask and correct aspect ratio
        plt.figure()
        x=sns.heatmap(similarity_matrix_plot, mask=mask, cmap=cmap, vmin=0.55, vmax=0.75,center=0.65,
                       xticklabels=similarity_matrix.columns,yticklabels=similarity_matrix.columns,
                       square=True, linewidths=.5, annot_kws={"size": 16}, cbar_kws={"shrink": .6})
        #ax[c//2, c%2].set_title('{}'.format(selected_chr))
        x.set_title('{}'.format(selected_chr))
        #outdir='/Users/jianingliu/Downloads/'
        plt.savefig(os.path.join(outdir, '.'.join((selected_chr, "similarity.pdf"))), bbox_inches="tight", pad_inches=0)

    print("Permutation analysis")
    '''permutate the order of NAMs'''
    '''saturation calculation'''
    saturation_stat_allchrs=pd.DataFrame()
    perm = 1000
    for selected_chr in interval_pangenome_allchrs['chr'].unique():
        #selected_chr='chr10'
        interval_pangenome_binary=interval_pangenome_allchrs[interval_pangenome_allchrs['chr']==selected_chr]
        saturation_stat=pd.DataFrame()
        saturation_stat[0]=range(1,len(NAMline)+1)
        for i in range(0,perm):
            #i=0
            ordered=shuffled(NAMline)
            counted_len=0
            counted_index=[]
            for j in range(len(ordered)):
                #j=0
                line=ordered[j]
                selected=interval_pangenome_binary[~interval_pangenome_binary.index.isin(counted_index)]
                counted=selected[selected[line]==1]
                counted_len += (counted['end']-counted['start']+1).sum()
                counted_index+=counted.index.tolist()
                saturation_stat.loc[j,i+1]= counted_len
        saturation_stat=pd.DataFrame(saturation_stat).astype('int64')
        saturation_stat['chr']=selected_chr
        saturation_stat_allchrs=saturation_stat_allchrs.append(saturation_stat)
    '''Write permutation summary'''
    #saturation_stat_file='/Users/jianingliu/Downloads/sat_stat.txt'
    pd.DataFrame(saturation_stat_allchrs).to_csv(saturation_stat_file,index=False,header=False,sep='\t')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Frequency counts of ref genomes')
    parser.add_argument('-i','--alignment_file', required=True, help='Syntenic pairwise alignment file (after chaining)')
    parser.add_argument('-g','--genome_file', required=True,
                        help='chromosome size of input genomes for pan-genome calculation')
    parser.add_argument('-r','--ref', required=True)
    parser.add_argument('-d','--outdir',required=True)
    parser.add_argument('-f','--frequency_file',required=True,
                        help='output file for frequency of all segments across all chrs')
    parser.add_argument('-s','--saturation_stat_file',required=True,
                        help='permutation result of 1000 times, count the coverage of ref with alignment results')

    # parser.add_argument('-r','--ref', default='B73',
    #                     help='Ref genome ID')
    # parser.add_argument('-q','--query',
    #                     help='Query genome ID')

    args = parser.parse_args()
    main(alignment_file=args.alignment_file,genome_file=args.genome_file,ref=args.ref,outdir=args.outdir,frequency_file=args.frequency_file,saturation_stat_file=args.saturation_stat_file)






#for i in range(len(NAMline)):
#i=0
#NAMline_query=NAMline[:i] + NAMline[i+1:]
# for j in range(n_elements):
#     NAMline_permutation.append(shuffled(NAMline))
   # NAMline_permutation += [[NAMline[i]] + list(x) for x in NAMline_query_permutation]

        # sns.heatmap(similarity_matrix,
        #         xticklabels=similarity_matrix.columns,
        #         yticklabels=similarity_matrix.columns, vmin=0.5, vmax=0.7,center=0.6)



# import scipy.cluster.hierarchy as shc
# plt.figure(figsize=(10, 7))
# plt.title("Dendrograms")
# dend = shc.dendrogram(shc.linkage(similarity_matrix, method='ward'),labels=similarity_matrix.columns)
#interval_pangenome_binary_pca=interval_pangenome_binary[NAMline].T
#from sklearn.decomposition import PCA
#pca = PCA(interval_pangenome_binary_pca)
