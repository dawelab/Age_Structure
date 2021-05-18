#/usr/bin/env python3
import numpy as np
import pandas as pd
import ntpath
import os
import sys
import subprocess

#paf_file="/scratch/jl03308/NAM_pancentromere/genome_alignment/chr1_shujun/P39_chr1.mapped-to.B73_chr1.sorted.paf"
#syn_bed_file="/scratch/jl03308/NAM_pancentromere/NAM_SV/chr1/B73_P39.aligned.bed"
#vcf_file="/scratch/jl03308/NAM_pancentromere/NAM_SV/chr1/snp_window/P39_chr1.mapped-to.B73_chr1-199591_206534.syn.vcf"
def main(ref,paf_file,syn_bed_file,start_coords,end_coords,qstart_coords,qend_coords,output):

    chrs=ntpath.basename(paf_file).split(".")[0].split("_")[1]
    ref=ntpath.basename(paf_file).split(".")[2].split("_")[0]
    query=ntpath.basename(paf_file).split(".")[0].split("_")[0]
    start_coords=int(start_coords)
    end_coords=int(end_coords)

    wdir=ntpath.dirname(paf_file)
    path = wdir
    os.chdir(path)
    paf = pd.read_csv(paf_file,delimiter='\t',header=None,error_bad_lines=False, warn_bad_lines=False)
    # syn_paf_intersect_index=[]
    # for i in syn_selected.index:
    #     print(i)
    #     syn_ref_start=syn_selected.loc[i,1]
    #     syn_ref_end=syn_selected.loc[i,2]
    #     syn_query_start=syn_selected.loc[i,4]
    #     syn_query_end=syn_selected.loc[i,5]
    #     for j in paf_noseq_selected.index:
    #         paf_noseq_selected[(paf_noseq_selected[7]==syn_ref_start)&(paf_noseq_selected[8]==syn_ref_end)&(paf_noseq_selected[2]==syn_query_start)&(paf_noseq_selected[3]==syn_query_end)].index
    #         paf_ref_start=paf_noseq_selected.loc[j,7]
    #         paf_ref_end=paf_noseq_selected.loc[j,8]
    #         paf_query_start=paf_noseq_selected.loc[j,2]
    #         paf_query_end=paf_noseq_selected.loc[j,3]
    #         if syn_ref_start==paf_ref_start and syn_ref_end==paf_ref_end and syn_query_start==paf_query_start and syn_query_end==paf_query_end:
    #             syn_paf_intersect_index.append(j)

    subprocess.call(['sh', '/home/jl03308/git/NAM_pancentromere/NAM_SV/paftools_sv_call_ref.sh', str(chrs), str(ref),paf_file])
    #os.system("sh /home/jl03308/git/NAM_pancentromere/NAM_SV/paftools_sv_call.sh " + str(chrs) + " " + '.'.join(ntpath.basename(paf_file).split(".")[:3]) + "-" + str(start_coords) + "_" + str(end_coords) + '.syn.paf')

    vcf_file='.'.join(ntpath.basename(paf_file).split(".")[:3]) + '.syn.vcf'
    aligned_len=paf[10].sum()
    #if os.path.exists(vcf_file):
    try:
        vcf = pd.read_csv(vcf_file,delimiter='\t',header=None,comment='#')
    except pd.errors.EmptyDataError:
        vcf = pd.DataFrame()
    if not vcf.empty:
        vcf_snp_index=[]
        for i in vcf.index:
            if len(vcf.loc[i,3])==1 and len(vcf.loc[i,4])==1 and vcf.loc[i,3]!=vcf.loc[i,4]:
                vcf_snp_index.append(i)
        vcf.loc[vcf_snp_index].to_csv('.'.join(ntpath.basename(paf_file).split(".")[:3]) + "-" + str(start_coords) + "_" + str(end_coords) + '.syn.snp.vcf', sep='\t', index=False,header=False)
        divergence=len(vcf_snp_index)/aligned_len/2/3.3 *100000000
        with open(output, "a") as div:
            div.write(ref + "\t" + str(start_coords) +  "\t" + str(end_coords) + "\t" + query + "\t" + str(qstart_coords) +  "\t" + str(qend_coords) + "\t" + str(len(vcf_snp_index)) + "\t" + str(aligned_len) + "\t" + str(round(divergence,1)) + "\n")
    else:
        with open(output, "a") as div:
            div.write(ref + "\t" + str(start_coords) +  "\t" + str(end_coords) + "\t" + query + "\t" + str(qstart_coords) +  "\t" + str(qend_coords) + "\t" + str(0)  + "\t" + str(aligned_len)+ "\t" + str(0) + "\n")

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='SNP dating with paftools between two haplotypes')
    parser.add_argument('-r','--ref', required=True)
    parser.add_argument('-i1','--paf_file', required=True,
                        help='pairwise alignment output of minimap2')
    parser.add_argument('-i2','--syn_bed_file', required=True,
                        help='syntenic aligned bed file')
    parser.add_argument('-s','--start_coords',required=True,
                        help='start coordinate of selected region')
    parser.add_argument('-e','--end_coords',required=True,
                        help='end coordinate of selected region')
    parser.add_argument('-qs','--qstart_coords',required=True,
                        help='start coordinate of selected region of query')
    parser.add_argument('-qe','--qend_coords',required=True,
                        help='end coordinate of selected region of query')
    parser.add_argument('-o','--output',required=True,
                        help='output file')

    args = parser.parse_args()
    main(ref=args.ref, paf_file=args.paf_file, syn_bed_file=args.syn_bed_file,start_coords=args.start_coords,end_coords=args.end_coords,qstart_coords=args.qstart_coords,qend_coords=args.qend_coords,output=args.output)
