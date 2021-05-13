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

''' Two rounds of chaining, first chaining takes mapping quality into consideration,
    second round does not and aim to further chain the >15kb blocks missed by the first round'''

def chaining(query,alignedlen,mapqual):
    n = len(query)
    # Declare the list (array) for LIS and
    # initialize LIS values for all indexes
    score = [0]*n
    score[0] = alignedlen[0] * math.log10(mapqual[0] + 0.001) - abs(query[0]-0)/100
    prev = [0]*n
    for i in range(0, n):
        prev[i] = i

    # Compute optimized LIS values in bottom up manner
    for i in range (1 , n):
        for j in range(0 , i):
            if query[i] > query[j] and score[i] < score[j] + alignedlen[i] * math.log10(mapqual[i] + 0.001) - abs(query[i]-query[j])/100:
                score[i] = score[j]+ alignedlen[i] * math.log10(mapqual[i] + 0.001) - abs(query[i]-query[j])/100
                #aligned segment length * adjusted mapping quality value - adjusted gap penalty
                prev[i] = j
                #print([i,score[i]])
    # traceback
    #find the index of largest value
    idx=np.argmax(score)
    # Create variables to store index
    syn_index=[idx]
    while idx != prev[idx] and idx != 0:
        idx = prev[idx]
        syn_index.append(idx)
    #len(syn_index)
    return syn_index[::-1]

def chaining_nomapq(query,alignedlen):
    n = len(query)
    # Declare the list (array) for LIS and
    # initialize LIS values for all indexes
    score = [0]*n

    score[0] = alignedlen[0] - abs(query[0]-0)/100
    prev = [0]*n
    for i in range(0, n):
        prev[i] = i

    # Compute optimized LIS values in bottom up manner
    for i in range (1 , n):
        for j in range(0 , i):
            if query[i] > query[j] and score[i] < score[j] + alignedlen[i] - abs(query[i]-query[j])/100:
                score[i] = score[j]+ alignedlen[i] - abs(query[i]-query[j])/100
                prev[i] = j
    # traceback
    #find the index of largest value
    # Create variables to store index
    idx=np.argmax(score)
    syn_index=[idx]
    while idx != prev[idx]:
        idx = prev[idx]
        syn_index.append(idx)
    #len(syn_index)
    return syn_index[::-1]


'''cluster function to group by interval defined by distance, similar to bedtools merge'''
def interval_cluster(data,distance):
    index=list(np.argsort(data))
    data.sort()
    groups = [[data[0]]]
    group_index=[[index[0]]]
    for i in range(1,len(data)):
        if abs(data[i] - groups[-1][-1]) <= int(distance):
            groups[-1].append(data[i])
            group_index[-1].append(index[i])
        else:
            groups.append([data[i]])
            group_index.append([index[i]])
    return groups, group_index

def segment_cluster(data):
    cluster=1
    data.loc[0,'cluster']=1
    for i in data.index[1:]:
        if data.loc[i,12] == data.loc[i-1,12]:
            data.loc[i,'cluster'] = int(cluster)
        else:
            cluster+=1
            data.loc[i,'cluster'] = int(cluster)
    return data

'''mark duplicated region in syntenic aligned area'''

def dup_extract(data,cutoff):
    #data=merged1
    #cutoff=10000
    #syn_cluster=merged_syntenic.loc[1727:1743,[2,3,4,7,8,9,10,12,13]].reset_index(drop=True)
    #syn_cluster=merged_syntenic.loc[0:30,[2,3,4,7,8,9,10,12,13]].reset_index(drop=True)
    dup=pd.DataFrame()
    for cluster in data['cluster'].unique():
        #cluster=757
        syn_cluster= data[data['cluster']==cluster].reset_index(drop=True)
        # syn_cluster.loc[:,[2,3,7,8,10,11,12,13]]
        if syn_cluster[12].unique()[0] =="syn" and syn_cluster.shape[0]>1:
            for i in syn_cluster.index[1:]:
                #i=10
                if syn_cluster.loc[i-1,8] - syn_cluster.loc[i,7] > cutoff and syn_cluster.loc[i-1,3] < syn_cluster.loc[i,2]:
                    #print(cluster)
                    #print(i)
                    #syn_cluster.loc[i]
                    dup1=pd.Series(syn_cluster.loc[i,[0,1]].tolist() + [syn_cluster.loc[i-1,3]-abs(syn_cluster.loc[i-1,8]-syn_cluster.loc[i,7]),syn_cluster.loc[i-1,3]] + syn_cluster.loc[i,[4,5,6]].tolist() + [syn_cluster.loc[i,7],syn_cluster.loc[i-1,8]] + syn_cluster.loc[i,[9,10,11]].tolist() + ['tandemdup_q',syn_cluster.loc[i,13],'cluster'])
                    dup2=pd.Series(syn_cluster.loc[i,[0,1]].tolist() + [syn_cluster.loc[i,2],syn_cluster.loc[i,2]+abs(syn_cluster.loc[i-1,8]-syn_cluster.loc[i,7])]+ syn_cluster.loc[i,[4,5,6]].tolist() + [syn_cluster.loc[i,7],syn_cluster.loc[i-1,8]] + syn_cluster.loc[i,[9,10,11]].tolist() + ['tandemdup_q',syn_cluster.loc[i,13],'cluster'])
                    dup=dup.append(dup1,ignore_index=True)
                    dup=dup.append(dup2,ignore_index=True)
                #     #if abs(syn_cluster.loc[i-1,8]-syn_cluster.loc[i,7]) < min(syn_cluster.loc[i-1,10],syn_cluster.loc[i,10]) * 0.5:
                if syn_cluster.loc[i-1,3] - syn_cluster.loc[i,2] > cutoff and syn_cluster.loc[i-1,8] < syn_cluster.loc[i,7]:
                    #print(cluster)
                    #print(i)
                    #syn_cluster.loc[i]
                    dup1=pd.Series(syn_cluster.loc[i,[0,1]].tolist() + [syn_cluster.loc[i,2],syn_cluster.loc[i-1,3]] + syn_cluster.loc[i,[4,5,6]].tolist() + [syn_cluster.loc[i-1,8]- abs(syn_cluster.loc[i-1,3]-syn_cluster.loc[i,2]), syn_cluster.loc[i-1,8]] + syn_cluster.loc[i,[9,10,11]].tolist() + ['tandemdup_r',syn_cluster.loc[i,13],'cluster'])
                    dup2=pd.Series(syn_cluster.loc[i,[0,1]].tolist() + [syn_cluster.loc[i,2],syn_cluster.loc[i-1,3]] + syn_cluster.loc[i,[4,5,6]].tolist() + [syn_cluster.loc[i,7],syn_cluster.loc[i,7] + abs(syn_cluster.loc[i-1,3]-syn_cluster.loc[i,2])] + syn_cluster.loc[i,[9,10,11]].tolist() + ['tandemdup_r',syn_cluster.loc[i,13],'cluster'])
                    dup=dup.append(dup1,ignore_index=True)
                    dup=dup.append(dup2,ignore_index=True)
    if dup.shape[0] >=1:
        dup[10]=abs(dup[2]-dup[3])
    #dup.loc[:,[2,3,7,8,10,11,12,13]]
        return dup
    else:
        return None



def main(paf_file,syntenic_alignment_out_file):

    '''read file and filter'''
    aligned_segments=pd.read_csv(paf_file, delimiter="\t",comment="#",header=None).sort_values([2,3],ascending=[True,False]).reset_index(drop=True).loc[:,0:11]
    '''remove complete embed'''
    embed_index=[]
    x=0
    for i in aligned_segments.index[1:]:
        if aligned_segments.loc[i,2] <= aligned_segments.loc[x,3] and aligned_segments.loc[i,3] <= aligned_segments.loc[x,3]: # complete embedment
            embed_index.append(i)
        else:
            x=i
    aligned_segments_noembed=aligned_segments.drop(embed_index).sort_values([7,8],ascending=[True,False]).reset_index(drop=True)
    aligned_segments=aligned_segments_noembed

    embed_index=[]
    x=0
    for i in aligned_segments.index[1:]:
        if aligned_segments.loc[i,7] <= aligned_segments.loc[x,8] and aligned_segments.loc[i,8] <= aligned_segments.loc[x,8]: # complete embedment
            embed_index.append(i)
        else:
            x=i
    aligned_segments_noembed=aligned_segments.drop(embed_index).sort_values([7,8],ascending=[True,False]).reset_index(drop=True)
    aligned_segments=aligned_segments_noembed
    #aligned_segments=aligned_segments.head(500)
    '''chaining'''
    '''divide into inv and noninv groups'''
    print("1. Performing the first round of chaining")
    noninv=aligned_segments[aligned_segments[4]=="+"]
    #synaligned_noninv.loc[16:28,[2,3,7,8,10,11,4]]
    ref=noninv[7].to_numpy()
    query=noninv[2].to_numpy()
    alignedlen=noninv[10].to_numpy()
    mapqual=noninv[11].to_numpy()

    syn_index=chaining(query,alignedlen,mapqual)
    #len(syn_index)
    #as the order was calculated by chaining, index needs to be projected onto selected index
    syn_noninv_index=noninv.index[syn_index]
    #syn_noninv_index[0:20]
    synaligned_noninv=noninv.loc[syn_noninv_index]
    #synaligned[synaligned[4]=="-"]
    #synaligned_noninv.loc[:,[2,3,7,8,4,10,11]]
    #noninv.loc[:,[2,3,7,8,4,10,11]]

    '''inversion'''
    '''steps different from syntenic segment chaining, as inversions are distributed in certain area of genome'''
    inv=aligned_segments[aligned_segments[4]=="-"]
    '''remove inv that overlap with synaligned_noninv dataframe'''
    overlap_index=[]
    for i in inv.index:
        inv_interval=pd.Interval(inv.loc[i,7], inv.loc[i,8])
        for j in synaligned_noninv.index:
            synaligned_noninv_interval=pd.Interval(synaligned_noninv.loc[j,7], synaligned_noninv.loc[j,8])
            if inv_interval.overlaps(synaligned_noninv_interval):
                if i not in overlap_index:
                        overlap_index.append(i)
    #len(overlap_index)
    inv=inv.drop(overlap_index)
    #inv.loc[(inv[2]>99000000)&(inv[2]<102000000)&(inv[7]>99000000)&(inv[7]<103000000),[2,3,7,8,10,11]]
    '''chain inversion'''
    print("Detecting inversions")
    #break inversions by index
    #real ancient inversions should be consisted of arrays of small conjacent inversions
    index_breaks=interval_cluster(inv.index.tolist(),2)[0]
    syninv_index_pass1=[]
    #size=[]
    for x in index_breaks:
        #if each inversion set consists of more than 2 segments, considered a true inversion
        if len(x) > 2:
            #print(x)
            #x=[12074, 12076, 12077, 12078, 12079, 12081]
            pre_group=aligned_segments.loc[x,2].values.tolist()
            #break inversion by interval distance
            div_group=interval_cluster(pre_group,500000)[1]
            div_group_index=[[aligned_segments.loc[x,2].index[item] for item in list] for list in div_group]
            #chain each inversion segments
            for y in div_group_index:
                if len(y) > 2:
                    #print(y)
                    #y=[12079, 12078, 12076]
                    #sort values by opposite order as positive orientation
                    inv_tmp=aligned_segments.loc[y].sort_values([7],ascending=[False])
                    query=inv_tmp[2].to_numpy()
                    alignedlen=inv_tmp[10].to_numpy()
                    mapqual=inv_tmp[11].to_numpy()
                    inv_pass1_index=inv_tmp.index[chaining(query,alignedlen,mapqual)].tolist()[::-1]
                    syninv_index_pass1.append(inv_pass1_index)
                    #if aligned_segments.loc[inv_pass1_index][10].sum()> 50000:
                        #print(aligned_segments.loc[inv_pass1_index,[2,3,7,8,10,11,4]])
                    #as the order was calculated by chaining, index needs to be projected onto inversion index
                    #size.append(aligned_segments.loc[inv_pass1_index][10].sum())
    #len(syninv_index_pass1)
    synaligned_inv=inv.loc[[item for sublist in syninv_index_pass1 for item in sublist]]
    for i in range(len(syninv_index_pass1)):
        synaligned_inv.loc[syninv_index_pass1[i],13] = int(i+1)

    #synaligned_inv.loc[:,[2,3,7,8,13]]
    #synaligned_inv.loc[(synaligned_inv[2]>99000000)&(synaligned_inv[2]<102000000)&(synaligned_inv[7]>99000000)&(synaligned_inv[7]<103000000),[2,3,7,8,10,11,13]]
    '''other large segments that were not chained as inversion or syntenic alignment'''
    '''could be translocation or inv/syn missed by the first step because of low mapq'''
    '''chain the remaining segments'''
    #syn_index[10:30]
    #len(syn_noninv_index.tolist() + [item for sublist in syninv_index_pass1 for item in sublist])
    other=aligned_segments.drop(syn_noninv_index.tolist() + [item for sublist in syninv_index_pass1 for item in sublist])
    '''noninv'''
    print("Detecting translocations")
    other_noninv=other[other[4]=="+"]
    index_breaks=interval_cluster(other_noninv.index.tolist(),2)[0]
    other_noninv_index_pass1=[]
    #size=[]
    for x in index_breaks:
        #if each inversion set consists of more than 2 segments, considered a true inversion
        if len(x) > 2:
            #print(x)
            #x=[13403, 13404, 13405, 13406, 13407]
            pre_group=aligned_segments.loc[x,2].values.tolist()
            #break inversion by interval distance
            div_group=interval_cluster(pre_group,500000)[1]
            div_group_index=[[aligned_segments.loc[x,2].index[item] for item in list] for list in div_group]
            #chain each inversion segments
            for y in div_group_index:
                if len(y) > 2:
                    #print(y)
                    #y=[9038, 9037, 9036]
                    #other.loc[y].loc[:,[2,3,7,8,10,11]]
                    #sort values by opposite order as positive orientation
                    other_noninv_tmp=other.loc[y].sort_values([7],ascending=[True])
                    query=other_noninv_tmp[2].to_numpy()
                    alignedlen=other_noninv_tmp[10].to_numpy()
                    mapqual=other_noninv_tmp[11].to_numpy()
                    other_noninv_tmp_index=other_noninv_tmp.index[chaining_nomapq(query,alignedlen)].tolist()[::-1]
                    other_noninv_index_pass1.append(other_noninv_tmp_index)
    other_noninv_syn=other.loc[[item for sublist in other_noninv_index_pass1 for item in sublist]]
    for i in range(len(other_noninv_index_pass1)):
        other_noninv_syn.loc[other_noninv_index_pass1[i],13] = int(i+1)
    '''inv'''
    print("Detecting inverted translocations")
    other_inv=other[other[4]=="-"]
    other_inv_index_pass1=[]
    index_breaks=interval_cluster(other_inv.index.tolist(),2)[0]
    #size=[]
    for x in index_breaks:
        #if each inversion set consists of more than 2 segments, considered a true inversion
        if len(x) > 2:
            #print(x)
            #x=[13403, 13404, 13405, 13406, 13407]
            pre_group=aligned_segments.loc[x,2].values.tolist()
            #break inversion by interval distance
            div_group=interval_cluster(pre_group,500000)[1]
            div_group_index=[[aligned_segments.loc[x,2].index[item] for item in list] for list in div_group]
            #chain each inversion segments
            for y in div_group_index:
                if len(y) > 2:
                    #print(y)
                    #y=2310, 2311, 2312]
                    #other.loc[y].loc[:,[2,3,7,8,10,11]]
                    #sort values by opposite order as positive orientation
                    other_inv_tmp=other.loc[y].sort_values([7],ascending=[False])
                    query=other_inv_tmp[2].to_numpy()
                    alignedlen=other_inv_tmp[10].to_numpy()
                    mapqual=other_inv_tmp[11].to_numpy()
                    other_inv_tmp_index=other_inv_tmp.index[chaining_nomapq(query,alignedlen)].tolist()[::-1]
                    other_inv_index_pass1.append(other_inv_tmp_index)
    other_inv_syn=other.loc[[item for sublist in other_inv_index_pass1 for item in sublist]]
    #other_inv_syn.loc[:,[2,3,7,8]]
    for i in range(len(other_inv_index_pass1)):
        other_inv_syn.loc[other_inv_index_pass1[i],13] = int(i+1) + len(syninv_index_pass1)

    '''sizecutoff'''
    other_nonsyn=other.drop([item for sublist in other_inv_index_pass1 for item in sublist] + [item for sublist in other_noninv_index_pass1 for item in sublist])
    cutoff=15000
    other_size_cutoff=other_nonsyn[other_nonsyn[10]>cutoff]

    '''concatenate syn, inv, and other'''
    '''label'''
    synaligned_noninv[12] = "syn"
    true_inv_size=abs(synaligned_inv.loc[synaligned_inv[10] > 15000,2] - synaligned_inv.loc[synaligned_inv[10] > 15000,7])
    up=true_inv_size.mean() + 2*true_inv_size.std()
    for i in synaligned_inv[13].unique():
        #i=10
        interval=abs(synaligned_inv.loc[synaligned_inv[13]==i,2] - synaligned_inv.loc[synaligned_inv[13]==i,7]).mean()
        if interval > min(50000000,max(up,10000000)):
            synaligned_inv.loc[synaligned_inv[13]==i,12] = "trans_inv"
        else:
            synaligned_inv.loc[synaligned_inv[13]==i,12] = "inv"
    #synaligned_inv[synaligned_inv[12] == "inv"].loc[:,[2,3,7,8,10,13]]
    #        synaligned_inv.loc[:,[2,3,7,8]]
    #other_inv_syn[other_inv_syn[12] == "trans_inv"].loc[:,[2,3,7,8]]
    other_inv_syn[12] = "trans_inv"
    other_noninv_syn[12] = "trans"
    other=other_size_cutoff.copy()
    other[12] = "other"

    '''merge all three dataframe'''
    merged1=synaligned_noninv.append(synaligned_inv).append(other_inv_syn).append(other_noninv_syn).append(other).sort_values([7,2],ascending=True).reset_index(drop=True)
    #merged1[merged1[12]=="other"]
   #merged1.loc[7600:7620,[2,3,7,8,10,11,12]]
    print("2. Performing the second round of chaining")

    '''second round of chaining'''
    '''chaining without mapping quality'''
    '''The aim of second-round chaining is to chain the segments missed by the first round,
    (likely due to low mapping quality, labeled as "other"), they'd be relabelled as inv or syn after this step'''
    '''noninv'''
    #merged1_noninv=merged1.loc[merged1[12].str.contains("syn")]
    merged1_noninv=merged1[(merged1[12]=="syn") | (merged1[12]=="other") &(merged1[4]=="+")]
    #merged1_noninv=merged1_noninv.loc[2955:2965,[2,3,7,8,10,11,12,4]]
    #synaligned_noninv.loc[16:28,[2,3,7,8,10,11,4]]
    query=merged1_noninv[2].to_numpy()
    alignedlen=merged1_noninv[10].to_numpy()
    syn_index=chaining_nomapq(query,alignedlen)
    merged1.loc[merged1_noninv.index[syn_index],12] = "syn"
    #merged1.loc[(merged1[4]=="-")&(merged1[12]=="syn"),12] = "inv"
    #merged1.loc[(merged1[12]=="other")&(merged1[4]=="+")]
    '''inv'''
    merged1_inv=merged1[(merged1[12].str.contains("inv")) | (merged1[12]=="other") &(merged1[4]=="-")]#.sort_values([7],ascending=[False])
    index_breaks=interval_cluster(merged1_inv.index.tolist(),1)[0]
    for x in index_breaks:
        if merged1_inv.loc[x,12].str.contains('other').any() and len(x) > 2:
            #print(x)
            #print(merged1_inv.loc[:,[2,3,7,8,10,11,4,12,13]])
    #         #x=[6646, 6647, 6648, 6649, 6650, 6651, 6652, 6653, 6654]
            #x=[3705, 3706, 3707, 3708, 3709, 3710, 3711, 3712, 3713, 3714]
            pre_group=merged1.loc[x,2].values.tolist()
            div_group=interval_cluster(pre_group,500000)[1]
            div_group_index=[[aligned_segments.loc[x,2].index[item] for item in list] for list in div_group]
            for y in div_group_index:
                if merged1_inv.loc[y,12].str.contains('other').any() and len(y) > 1:
                    #print(y)
                    #print(merged1.loc[y,[2,3,7,8,10,11,4,12,13]])
    #               #y=[6654, 6653, 6652, 6648, 6647, 6646]
                    inv_tmp=merged1.loc[y].sort_values([7],ascending=[False])
                    query=inv_tmp[2].to_numpy()
                    alignedlen=inv_tmp[10].to_numpy()
                    inv_pass1_index=inv_tmp.index[chaining_nomapq(query,alignedlen)].tolist()[::-1]
                    if merged1.loc[inv_pass1_index,13].dropna().unique().shape[0]>0:
                        inv_index=merged1.loc[inv_pass1_index,13].dropna().unique()[0]
                    else:
                        inv_index=merged1[13].max() +1
                    merged1.loc[inv_pass1_index,12] = 'inv'
                    merged1.loc[inv_pass1_index,13] = inv_index
    #merged2.loc[(merged2[2]>99000000)&(merged2[2]<99123778)&(merged2[7]>99000000)&(merged2[7]<103000000),[2,3,7,8,10,11,13]]
    #merged2.loc[(merged2[13]==122)|(merged2[13]==123),[2,3,7,8,4,10,12,13]]

    '''segment all syntenic blocks into clusters'''
    merged1=segment_cluster(merged1)
    cutoff=10000
    print("Detecting tandem duplications")
    '''identify tandem duplications (>10Kb) in each syntenic cluster'''
    merged1_dup=dup_extract(merged1,cutoff)
    #merged1_dup.loc[:,[2,3,7,8,10,11]]
    merged2=merged1.append(merged1_dup).sort_values([7,2],ascending=True).reset_index(drop=True)
    #merged2[merged2[12].str.contains('trans')].loc[:,[2,3,7,8,13]]
    #merged2.loc[merged2[12].str.contains('inv'),[2,3,7,8,13]]
    #merged2[merged2['size'] ==414486].loc[:,[2,3,7,8,'size']]
    print("Calculating total size of anchors in each inversion and translocation")
    '''label the inversions and translocations by size, could be used for downstream filtering'''
    #inv_sizecutoff=20000
    #filtered_index=[]
    #all_inv_sizecutoff=pd.DataFrame()
    for i in merged2.loc[merged2[12].str.contains('inv'),13].dropna().unique():
        #i=123
        merged2.loc[(merged2[12].str.contains('inv'))&(merged2[13]==i),'size'] = merged2.loc[(merged2[12].str.contains('inv'))&(merged2[13]==i),10].sum()
    for i in merged2.loc[merged2[12]=="trans",13].dropna().unique():
        merged2.loc[(merged2[12]=="trans")&(merged2[13]==i),'size'] = merged2.loc[(merged2[12]=="trans")&(merged2[13]==i),10].sum()

    merged2=segment_cluster(merged2)
    #'''adjust overlapping coordinates in ref in each syntenic cluster'''
    #merged2_coordadjusted, pre_unalignment, trans =coords_adjust(merged2)

    print("Writing aligned blocks (synteny, inversion, translocation)")

    '''output aligned blocks'''
    ref=aligned_segments[5].unique()[0].split('_')[0]
    query=aligned_segments[0].unique()[0].split('_')[0]
    #merged2_sum.loc[0]
    merged2[14]=merged2[9]/merged2[10]*100
    merged2[15]=ref
    merged2[16]=query
    merged2_sum=merged2.loc[:,[15,7,8,16,2,3,4,14,12,9,10,13,'size']]
    merged2_sum['size']=merged2_sum['size'].fillna(0)
    merged2_sum[13]=merged2_sum[13].fillna(0)
    merged2_sum.loc[:,[7,8,2,3,14,9,10,13,'size']]=merged2_sum.loc[:,[7,8,2,3,14,9,10,13,'size']].astype(int)
    merged2_sum.loc[:,[15,7,8,16,2,3,4,10,14,12,13,'size']].to_csv(syntenic_alignment_out_file,sep="\t",index=False,header=False)
    #merged2_sum[merged2_sum[12]=="inv"]

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='extract syntenic alignment and unaligned segments from paf file')
    parser.add_argument('-i','--paf_file', required=True, help='paf pairwise alignment file (output of minimap2)')
    parser.add_argument('-o','--syntenic_alignment_out_file', required=True,
                        help='syntenic aligned segments between two input genomes')

    args = parser.parse_args()
    main(paf_file=args.paf_file,syntenic_alignment_out_file=args.syntenic_alignment_out_file)
