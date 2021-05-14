#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from itertools import permutations, islice
import random
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from interval import Interval, IntervalSet

# flint_origin=['HP301','P39','Il14H']
# temperate_origin=['B97','Ky21','M162W','Ms71','Oh43','Oh7b']
# mixed_origin = ['M37W','Mo18W','Tx303']
# tropical_origin = ['CML52','CML69','CML103','CML228','CML247','CML277','CML322','CML333','Ki3','Ki11','NC350','NC358','Tzi8']
# NAMline_origin=flint_origin+temperate_origin+mixed_origin+tropical_origin


'''merge overlaps , similar to bedtools merge'''
def merge(df):
    #df=alignment_added_genome
    df=df.sort_values([1],ascending=[True]).reset_index(drop=True)
    size=[df.shape[0]]
    i=0
    x=0
    while x < size[i-1]:
        df["group"]=(df[1]>df[2].shift()).cumsum()
        df=df.groupby("group").agg({1:"min", 2: "max"})[[1,2]]
        df.index.name = None
        size.append(df.shape[0])
        x=df.shape[0]
        i+=1
    size=(df[2]-df[1] + 1).sum()
    return df, size


'''permutate the order of NAMs'''
def shuffled(array):
    random.shuffle(array)
    return list(array)

#NAMline=['B73','B97','CML228','CML247','CML277']
def main(alignment_file,output,permutation_time):
    #NAMline=['B73','B97','CML103','CML228','CML247','CML277','CML322','CML333','CML52','CML69','HP301','IL14H','Ki11','Ki3','Ky21','M162W','M37W','Mo18W','MS71','NC350','NC358','Oh43','Oh7b','P39','Tx303','Tzi8']
    '''read file'''
    #chrs='chr10'
    #alignment_file='/Users/jianingliu/Downloads/test.bed'
    #output='/Users/jianingliu/Downloads/sat.bed'
    aligned_segments=pd.read_csv(alignment_file, delimiter="\t",comment="#",header=None).sort_values([0,3,1]).reset_index(drop=True)
    querylines=aligned_segments[3].unique()
    pangenome=pd.DataFrame()
    pangenome[0]=range(1,len(querylines)+1)
    '''shuffle query genomes and calculate total alignment space'''
    #permutation_time=3
    permutation_time=int(permutation_time)
    for i in range(0,permutation_time):
        querylines=shuffled(querylines)
        shuf=[]
        for i in range(1,len(querylines)+1):
            #print(i)
            x=querylines[0:i]
            aligned=aligned_segments[aligned_segments[3].isin(x)].sort_values([0,1,2]).reset_index(drop=True)
            size=merge(aligned)[1]
            shuf.append(size)
        pangenome=pd.concat([pangenome, pd.Series(shuf)], axis=1)

    pangenome.to_csv(output,sep="\t",index=False,header=False)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='calculate pan-genome space based on nx(n-1)/2 pairwise-genome alignment')
    parser.add_argument('-i','--alignment_file', required=True, help='Syntenic pairwise alignment file (after chaining)')
    parser.add_argument('-o','--output',required=True)
    parser.add_argument('-n','--permutation_time',required=True)

    # parser.add_argument('-r','--ref', default='B73',
    #                     help='Ref genome ID')
    # parser.add_argument('-q','--query',
    #                     help='Query genome ID')

    args = parser.parse_args()
    main(alignment_file=args.alignment_file,output=args.output,permutation_time=args.permutation_time)
