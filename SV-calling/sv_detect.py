#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import math
import itertools


'''extract unaligned regions'''

'''1) cluster sytenic blocks by breakpoints (inv, translocation)'''
def segment_cluster(data):
    cluster=1
    data.loc[0,'cluster']=1
    for i in data.index[1:]:
        if data.loc[i,9] == data.loc[i-1,9]:
            data.loc[i,'cluster'] = int(cluster)
        else:
            cluster+=1
            data.loc[i,'cluster'] = int(cluster)
    return data


'''2) as there are overlaps in ref and query coords, adjust coords  according to ref
for example:
1827  119523529  122628577  +  chr5  226353449  119954051  123059106  3105387
1828  122625251  122861321  +  chr5  226353449  123048841  123284852   236110
adjusted according to column 7, the result is:
1827  119523529  122628577  +  chr5  226353449  119954051  123059106  3105387
1828  122635561  122861321  +  chr5  226353449  123059106  123284852   236110
    '''
#merged_syntenic.shape
def coords_adjust(data):
    #data=syn_aligned_segmented
    merged_syntenic_adjusted=pd.DataFrame()
    for cluster in data['cluster'].unique():
        #cluster=11
        syn_cluster= data[data['cluster']==cluster].reset_index(drop=True)
        # data.loc[3780:3790,[1,2,4,5,9]]
        if syn_cluster[9].unique()[0] =="syn" and syn_cluster.shape[0]>1:
            for i in syn_cluster.index[1:]:
                #i=12
                if syn_cluster.loc[i-1,2] > syn_cluster.loc[i,1]:
                    if syn_cluster.loc[i-1,5] > syn_cluster.loc[i,4]:
                        #print(syn_cluster.loc[i])
                        #syn_cluster.loc[i, 'type'] = 'repeat'
                        syn_cluster.loc[i-1, 'type'] = 'repeat'
                    else:
                        if syn_cluster.loc[i,7] > syn_cluster.loc[i-1,7]: # adjust the coordinates of segments of the smaller size
                            syn_cluster.loc[i-1,5] = syn_cluster.loc[i-1,5] - abs(syn_cluster.loc[i,1] - syn_cluster.loc[i-1,2])
                            syn_cluster.loc[i-1,2] =syn_cluster.loc[i,1]
                            syn_cluster.loc[i-1,7] = max(abs(syn_cluster.loc[i-1,2] - syn_cluster.loc[i-1,1]), abs(syn_cluster.loc[i-1,5] - syn_cluster.loc[i-1,4]))
                        else:
                            syn_cluster.loc[i,4] = syn_cluster.loc[i,4] + abs(syn_cluster.loc[i-1,2] - syn_cluster.loc[i,1])
                            syn_cluster.loc[i,1] = syn_cluster.loc[i-1,2]
                            syn_cluster.loc[i,7] = max(abs(syn_cluster.loc[i,2] - syn_cluster.loc[i,1]), abs(syn_cluster.loc[i,5] - syn_cluster.loc[i,4]))
                if syn_cluster.loc[i-1,5] > syn_cluster.loc[i,4] and syn_cluster.loc[i-1,2] <= syn_cluster.loc[i,1]:
                        if syn_cluster.loc[i,7] > syn_cluster.loc[i-1,7]: # adjust the coordinates of segments of the smaller size
                            syn_cluster.loc[i-1,2] = syn_cluster.loc[i-1,2] - abs(syn_cluster.loc[i,4] - syn_cluster.loc[i-1,5])
                            syn_cluster.loc[i-1,5] =syn_cluster.loc[i,4]
                            syn_cluster.loc[i-1,7] = max(abs(syn_cluster.loc[i-1,2] - syn_cluster.loc[i-1,1]), abs(syn_cluster.loc[i-1,5] - syn_cluster.loc[i-1,4]))
                        else:
                            syn_cluster.loc[i,1] = syn_cluster.loc[i,1] + abs(syn_cluster.loc[i-1,5] - syn_cluster.loc[i,4])
                            syn_cluster.loc[i,4] = syn_cluster.loc[i-1,5]
                            syn_cluster.loc[i,7] = max(abs(syn_cluster.loc[i,2] - syn_cluster.loc[i,1]), abs(syn_cluster.loc[i,5] - syn_cluster.loc[i,4]))

            merged_syntenic_adjusted=merged_syntenic_adjusted.append(syn_cluster,ignore_index=True)
#merged_syntenic_adjusted[merged_syntenic_adjusted['type']=='repeat'][[1,2,4,5]]
        elif syn_cluster[9].unique()[0] =="inv" and syn_cluster.shape[0]>1:
            for i in syn_cluster.index[1:]:
                #i=12
                if syn_cluster.loc[i-1,2] > syn_cluster.loc[i,1]:
                    if syn_cluster.loc[i,5] > syn_cluster.loc[i-1,4]:
                        syn_cluster.loc[i, 'type'] = 'repeat'
                        #syn_cluster.loc[i-1, 'type'] = 'repeat'
                    else:
                        if syn_cluster.loc[i,7] > syn_cluster.loc[i-1,7]: # adjust the coordinates of segments of the smaller size
                            syn_cluster.loc[i-1,4] = syn_cluster.loc[i-1,4] + abs(syn_cluster.loc[i,1] - syn_cluster.loc[i-1,2])
                            syn_cluster.loc[i-1,2] =syn_cluster.loc[i,1]
                            syn_cluster.loc[i-1,7] = max(abs(syn_cluster.loc[i-1,2] - syn_cluster.loc[i-1,1]), abs(syn_cluster.loc[i-1,5] - syn_cluster.loc[i-1,4]))
                        else:
                            syn_cluster.loc[i,5] = syn_cluster.loc[i,5] - abs(syn_cluster.loc[i-1,2] - syn_cluster.loc[i,1])
                            syn_cluster.loc[i,1] = syn_cluster.loc[i-1,2]
                            syn_cluster.loc[i,7] = max(abs(syn_cluster.loc[i,2] - syn_cluster.loc[i,1]), abs(syn_cluster.loc[i,5] - syn_cluster.loc[i,4]))
            merged_syntenic_adjusted=merged_syntenic_adjusted.append(syn_cluster,ignore_index=True)

    merged_syntenic_adjusted=merged_syntenic_adjusted.sort_values([1,2],ascending=True).reset_index(drop=True)
    #data[data['type']=='complicated']
    #data.loc[7040:7060,[2,3,7,8,10,11,12,'type']]
    return merged_syntenic_adjusted


#merged_syntenic[merged_syntenic[12]=="inv"].loc[1727:1747,[2,3,4,5,6,7,8,10]]
#merged_syntenic[merged_syntenic[2]>117075421].loc[:,[2,3,4,5,6,7,8,10]].head(20)
'''3) assign start and end unaligned segements'''
def unalignment_cal(data,ref_gsize,query_gsize):
    #data=merged2_coordadjusted
    #data.loc[:,[2,3,7,8,10,12]]
    #data.loc[data.index[-1],[2,3,7,8,10,12]]
    #unaligned_segments.loc[:,[2,3,7,8,10,12]]
    unaligned_segments=pd.DataFrame()
    '''start'''
    if data.loc[data.index[0],9] == "syn":
        unaligned_segments=unaligned_segments.append(pd.Series([1,data.loc[data.index[0],1],1,data.loc[data.index[0],4],"syn"]),ignore_index=True)
    elif data.loc[data.index[0],9] == "inv":
        unaligned_segments=unaligned_segments.append(pd.Series([1,data.loc[data.index[0],1],1,data[4].min(),"inv"]),ignore_index=True)
    '''end'''
    if data.loc[data.index[-1],9] == "syn":
        unaligned_segments=unaligned_segments.append(pd.Series([data.loc[data.index[-1],2],ref_gsize,data.loc[data.index[-1],5],query_gsize,"syn"]),ignore_index=True)
    elif data.loc[data.index[-1],9] == "inv":
        unaligned_segments=unaligned_segments.append(pd.Series([data.loc[data.index[-1],2],ref_gsize,data[5].max(),query_gsize,"inv"]),ignore_index=True)

    '''3) calculate unaligned regions in each cluster '''
    for cluster in data['cluster'].unique():
        #cluster=1
        syn_cluster= data[data['cluster']==cluster].reset_index(drop=True)
        #syn_cluster.loc[:,[2,3,7,8,10]]
        if syn_cluster[9].unique()[0] =="syn" and syn_cluster.shape[0]>1:
            #print(cluster)
            for i in syn_cluster.index[1:]:
                unaligned=[]
                unaligned.append(syn_cluster.loc[i-1,2])
                unaligned.append(syn_cluster.loc[i,1])
                unaligned.append(syn_cluster.loc[i-1,5])
                unaligned.append(syn_cluster.loc[i,4])
                unaligned.append("syn")
                unaligned_segments=unaligned_segments.append(pd.Series(unaligned),ignore_index=True)
        if syn_cluster[9].unique()[0]=="inv" and syn_cluster.shape[0]>1:
            for i in syn_cluster.index[1:]:
                unaligned=[]
                unaligned.append(syn_cluster.loc[i-1,2])
                unaligned.append(syn_cluster.loc[i,1])
                unaligned.append(syn_cluster.loc[syn_cluster.index[i],5])
                unaligned.append(syn_cluster.loc[syn_cluster.index[i-1],4])
                unaligned.append("inv")
                unaligned_segments=unaligned_segments.append(pd.Series(unaligned),ignore_index=True)
    return unaligned_segments.sort_values([0],ascending=True).reset_index(drop=True)


def main(syntenic_alignment_file,selectedsyn_file,ref_gsize,query_gsize,sv_sum_file,unalignment_out_file,fig_file):
    #syntenic_alignment_file='/Users/jianingliu/Downloads/B73_P39.aligned.bed'
    #query_gsize=300441852
    #ref_gsize=308452471
    '''read file and filter'''
    query_gsize=int(query_gsize)
    ref_gsize=int(ref_gsize)
    alignedsegments=pd.read_csv(syntenic_alignment_file, delimiter="\t",comment="#",header=None)
    #alignedsegments.head()
    #alignedsegments.loc[range(1490,1496)][[0,1,2,3,4,5,9,11]]
    #alignedsegments[alignedsegments[11]>=0][[9,11]].drop_duplicates().boxplot(by=9)
    #alignedsegments[(alignedsegments[9]=='tandemdup_r')][[0,1,2,3,4,5]]
    #alignedsegments[(alignedsegments[9].str.contains('trans'))&(alignedsegments[11]>=50000)][[0,1,2,3,4,5,9]]
    syn_aligned=alignedsegments[(alignedsegments[9]=='syn') | ((alignedsegments[9]=='inv')&(alignedsegments[11]>=20000)) | ((alignedsegments[9].str.contains('trans'))&(alignedsegments[11]>=50000))| ((alignedsegments[9].str.contains('tandem')))].sort_values([1,2],ascending=True).reset_index(drop=True)
    syn_aligned_segmented=segment_cluster(syn_aligned)
    '''summarize each cluster'''
    sv_sum=pd.DataFrame()
    for x in syn_aligned_segmented['cluster'].unique():
        #x=2
        sv=syn_aligned_segmented[syn_aligned_segmented['cluster']==x]
        svtype=[]
        svtype.append(sv[0].unique()[0])
        svtype.append(min(sv[1].tolist()+sv[2].tolist()))
        svtype.append(max(sv[1].tolist()+sv[2].tolist()))
        svtype.append(sv[3].unique()[0])
        svtype.append(min(sv[4].tolist()+sv[5].tolist()))
        svtype.append(max(sv[4].tolist()+sv[5].tolist()))
        svtype.append(sv[9].unique()[0])
        sv_sum=sv_sum.append(pd.Series(svtype),ignore_index=True)
    sv_sum[7]=abs(sv_sum[2]-sv_sum[1])
    sv_sum[8]=abs(sv_sum[5]-sv_sum[4])
    #sv_sum.loc[:,[1,2,4,5,6]]
    '''identfy unaligned segments between syn and rearranged segments'''
    unaligned_interval=pd.DataFrame()
        #print(cluster)
    for i in sv_sum.index[1:]:
        unaligned=[]
        unaligned.append(sv_sum.loc[i-1,2])
        unaligned.append(sv_sum.loc[i,1])
        unaligned.append(sv_sum.loc[i-1,5])
        unaligned.append(sv_sum.loc[i,4])
        unaligned.append("syn")
        unaligned_interval=unaligned_interval.append(pd.Series(unaligned),ignore_index=True)
    unaligned_interval=unaligned_interval[(unaligned_interval[1] - unaligned_interval[0] > -100) & (unaligned_interval[3] - unaligned_interval[2] > -100)]
    for i in unaligned_interval.index:
        if unaligned_interval.loc[i,1] < unaligned_interval.loc[i,0]:
            unaligned_interval.loc[i,1] = unaligned_interval.loc[i,0]
        if unaligned_interval.loc[i,3] < unaligned_interval.loc[i,2]:
            unaligned_interval.loc[i,3] = unaligned_interval.loc[i,2]

    '''adjust overlapping coordinates in ref in each syntenic cluster'''
    merged2_coordadjusted=coords_adjust(syn_aligned_segmented)
    #merged2_coordadjusted=merged2_coordadjusted[merged2_coordadjusted['type']!='repeat']
    merged2_coordadjusted=merged2_coordadjusted[(merged2_coordadjusted[2]>merged2_coordadjusted[1])&(merged2_coordadjusted[5]>merged2_coordadjusted[4])]
    #merged2_coordadjusted.loc[3770:3777,[1,2,4,5,6]]
    print("3. Extracting unaligned segments")
    '''extract unaligned segments within each syntenic cluster'''
    unaligned=unalignment_cal(merged2_coordadjusted,ref_gsize,query_gsize)
    #unaligned[unaligned[1] < unaligned[0]]
    unaligned=unaligned.append(unaligned_interval).sort_values([0],ascending=True).reset_index(drop=True)
    print("Filtering unaligned blocks")
    '''filter unaligned based on overlap with trans or inv, syn and size'''
    #filter_index=[]
    # trans_q=pd.arrays.IntervalArray.from_tuples([tuple(x) for x in trans[[2,3]].values])

    # for j in unaligned_filter.index:
    #     syn_q=pd.Interval(unaligned_filter.loc[j,0], unaligned_filter.loc[j,1])
    #     if any(trans_q.overlaps(syn_q)) and j not in filter_index:
    #         filter_index.append(j)
    # len(filter_index)
    filter_index=[]
    #unaligned_filter.loc[4100:4120,[0,1,2,3]]
    '''remove unaligned intervals that have alignments in ref and query (10kb alignment size as cutoff)'''
    #aligned_r=pd.arrays.IntervalArray.from_tuples([tuple(x) for x in syn_aligned[[1,2]].values])
   # aligned_q=pd.arrays.IntervalArray.from_tuples([tuple(x) for x in syn_aligned[[4,5]].values])
    unaligned_filter=unaligned[(unaligned[1]>=unaligned[0])&(unaligned[3]>=unaligned[2])]
    for j in unaligned_filter.index:
        #j=3043
        #unaligned_filter.loc[4103]
        syn_aligned_ref_check=syn_aligned[(syn_aligned[1]>=unaligned_filter.loc[j,0])&(syn_aligned[2]<=unaligned_filter.loc[j,1])]
        syn_aligned_query_check=syn_aligned[(syn_aligned[4]>=unaligned_filter.loc[j,2])&(syn_aligned[5]<=unaligned_filter.loc[j,3])]
        #syn_aligned_query_check.loc[4115:4120,[0,1,2,3,4,5,9]]
        if syn_aligned_ref_check[7].sum() > 10000 and j not in filter_index:
            filter_index.append(j)
        elif syn_aligned_query_check[7].sum() > 10000 and j not in filter_index:
            filter_index.append(j)

    #len(filter_index)
    unaligned_segments=unaligned_filter.drop(filter_index).sort_values([0],ascending=True).reset_index(drop=True)

    unaligned_segments[5]=unaligned_segments[1]-unaligned_segments[0]   # query unaligned len
    unaligned_segments[6]=unaligned_segments[3]-unaligned_segments[2]   # ref unaligned len
    unaligned_segments[7]=abs(unaligned_segments[5]-unaligned_segments[6]) # ref-query unaligned len difference
    #unaligned_segments[unaligned_segments[7]>5000000]
    #unaligned_segments=unaligned_segments[unaligned_segments[7]<5000000]
    unaligned_segments.loc[(unaligned_segments[6]==0)&(unaligned_segments[5]>0),8] = "ins_specifc"
    unaligned_segments.loc[(unaligned_segments[6]==0)&(unaligned_segments[5]<0),8] = "del_specifc"
    unaligned_segments.loc[(unaligned_segments[5]<=0)&(unaligned_segments[6]>0),8] = "del_specifc"
    unaligned_segments.loc[(unaligned_segments[5]>0)&(unaligned_segments[6]>0)&(unaligned_segments[6]>unaligned_segments[5]),8] = "del_unalignment"
    unaligned_segments.loc[(unaligned_segments[5]>0)&(unaligned_segments[6]>0)&(unaligned_segments[6]<unaligned_segments[5]),8] = "ins_unalignment"
    unaligned_segments.loc[(unaligned_segments[5]>0)&(unaligned_segments[6]>0)&(unaligned_segments[6]==unaligned_segments[5]),8] = "unalignment"
    unaligned_segments[9]=alignedsegments[0].unique()[0]
    unaligned_segments[10]=alignedsegments[3].unique()[0]
    unaligned_segments_sum=unaligned_segments.loc[:,[9,0,1,10,2,3,4,8,7,5,6]]

    segments_size=pd.DataFrame()
    segments_size['unaligned_ref']=abs(unaligned_segments_sum[1]-unaligned_segments_sum[0])
    segments_size['unaligned_query']=abs(unaligned_segments_sum[3]-unaligned_segments_sum[2])
    segments_size['aligned_ref']=abs(syn_aligned[2]-syn_aligned[1])
    segments_size['aligned_query']=abs(syn_aligned[5]-syn_aligned[4])
    size_plot=segments_size.boxplot()
    #size_plot.set_ylabel("size (bp)")
    fig=size_plot.get_figure()
    fig.savefig(fig_file)

    print("Writing unaligned blocks")

    '''output unaligned blocks used for pairwise realignment characterization'''

    #size_diff_cutoff=1000
    #realignment_segments=unaligned_segments_sum[unaligned_segments_sum[5]>size_diff_cutoff]
    unaligned_segments_sum.loc[:,[0,1,2,3,7,5,6]]=unaligned_segments_sum.loc[:,[0,1,2,3,7,5,6]].astype(int)
    unaligned_segments_sum.loc[:,[9,0,1,10,2,3,4,8,7]].to_csv(unalignment_out_file,sep="\t",index=False,header=False)

    sv_sum.loc[:,[1,2,4,5,7,8]]=sv_sum.loc[:,[1,2,4,5,7,8]].astype(int)
    sv_sum.to_csv(sv_sum_file,sep="\t",index=False,header=False)
    syn_aligned.to_csv(selectedsyn_file,sep="\t",index=False,header=False)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='extract syntenic alignment and unaligned segments from paf file')
    parser.add_argument('-i','--syntenic_alignment_file', required=True,
                        help='syntenic aligned segments between two input genomes')
    parser.add_argument('-a','--selectedsyn_file', required=True,
                        help='syntenic aligned segments selected (no translocation)')
    parser.add_argument('-gs1','--ref_gsize', required=True,
                        help='chr size of reference')
    parser.add_argument('-gs2','--query_gsize', required=True,
                        help='chr size of query')
    parser.add_argument('-s','--sv_sum_file',required=True,
                        help='alignment and rearrangement summary')
    parser.add_argument('-o','--unalignment_out_file',required=True,
                        help='syntenic unaligned segments between two input genomes')
    parser.add_argument('-f','--fig_file',required=True,
                        help='plot of syntenic aligned/unaligned segments between two input genomes')

    args = parser.parse_args()
    main(syntenic_alignment_file=args.syntenic_alignment_file,selectedsyn_file=args.selectedsyn_file,ref_gsize=args.ref_gsize,query_gsize=args.query_gsize,sv_sum_file=args.sv_sum_file,unalignment_out_file=args.unalignment_out_file,fig_file=args.fig_file)
