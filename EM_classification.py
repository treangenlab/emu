#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 16:17:47 2020

@author: kcurry
"""

import os
import argparse
import pathlib
import subprocess


import math
import pysam
#import pyfastx
import numpy as np
import pandas as pd
from collections import Counter 
from flatten_dict import unflatten, flatten 

from Bio import SeqIO
from ete3 import NCBITaxa


def get_align_stats(cigar_stats):
    '''Convert list of CIGAR stats containing NM to a dict containing mismatches
        
        cigar_stats: array length 11 containing CIGAR stats in order cigar_headers
        (default order of pysam .get_cigar_stats() func)    
        return: dict of cigar stats
    '''
    cigar_headers = ['M','I','D','N','S','H','P','=','X','B','NM']
    char_align_stats = dict(zip(cigar_headers, cigar_stats))
    nm = char_align_stats.pop('NM')
    char_align_stats['X'] = nm - char_align_stats['I'] - char_align_stats['D']
    return char_align_stats

def get_char_align_probabilites(bwa_sam):
    '''P(match), P(mismatch), P(insertion), P(deletion), P(softclipping), P(hardclipping) 
        by counting how often the corresponding operations occur in the primary alignments
        and by normalizing over the total sum of operations
    
        bwa_sam: path to sam file bwa output
        return: dict of likelihood for each cigar output with prob > 0
    '''
    # initialize character alignment probabilites
    samfile = pysam.AlignmentFile(bwa_sam)
    char_align_stats_primary = {}
    for read in samfile.fetch():
        if not read.is_secondary:
            read_align_stats = get_align_stats(read.get_cigar_stats()[0])
            char_align_stats_primary = Counter(char_align_stats_primary) + Counter(read_align_stats)
    n_char = sum(char_align_stats_primary.values())
    return {char:val/n_char for char, val in char_align_stats_primary.items()}

def compute_log_L_rgs(p_char, cigar_stats):
    ''' log(L(r|s)) = log(P(match)) ×n_matches + log(P(mismatches)) ×n_mismatches ...
    
        p_dict: dict of CIGAR likelihoods
        cigar_stats: dict of cigar stats to compute
        return: float log(L(r|s)) for sequences r and s in cigar_stats alignment
    '''
    value = 0
    for char in ['M','X','I','D','S','H']:
        if char not in p_char and cigar_stats[char] > 0:
            return None
        elif char in p_char:
            value += math.log(p_char[char])*cigar_stats[char]    
    return value

def log_L_rgs_dict(bwa_sam, p_char):
    '''dict containing log(L(r|s)) for all pairwise alignments in bwa output
    
        bwa_sam: path to sam file bwa output
        p_char: dict of likelihood for each cigar output with prob > 0
        return: dict[(r,s)]=log(L(r|s))
    '''
    # calculate log(L(r|s)) for all alignments
    log_L_rgs = {}
    samfile = pysam.AlignmentFile(bwa_sam)
    for alignment in samfile.fetch():
        if alignment.reference_name != None:
            val = compute_log_L_rgs(p_char, get_align_stats(alignment.get_cigar_stats()[0]))
            if ((alignment.query_name, alignment.reference_name) not in log_L_rgs or log_L_rgs[(alignment.query_name, alignment.reference_name)] < val) and val != None:
                log_L_rgs[(alignment.query_name, alignment.reference_name)] = val
    return log_L_rgs
            

def EM_iterations(log_L_rgs, db_ids):
    '''Expectation maximization algorithm for alignments in log_L_rgs dict
        
        log_L_rgs: dict[(r,s)]=log(L(r|s))
        db_ids: array of ids in database
        n_reads: int number of reads in input sequence (|R|)
        return: composition vector f[s]=relative abundance of s from sequence database 
    '''
    n_db = len(db_ids)
    f = dict.fromkeys(db_ids, 1/n_db)
    
    total_log_likelihood = -math.inf
    while (True):
        # log_L_rns = log(L(r|s))+log(f(s)) for each alignment
        log_L_rns = {}
        for (r,s) in log_L_rgs:
            if s in f and f[s] != 0: 
                log_L_rns[(r,s)] = (log_L_rgs[(r,s)]+math.log(f[s]))
        log_c = {r:-max(smap.values()) for r,smap in unflatten(log_L_rns).items()}
        L_rns_c = {(r,s):(math.e**(v+log_c[r])) for (r,s),v in log_L_rns.items()}
        L_r_c = {r:sum(s.values()) for r,s in unflatten(L_rns_c).items()}
    
        # P(s|r), likelihood read r emanates db seq s
        L_sgr = unflatten({(s,r):v/L_r_c[r] for (r,s),v in L_rns_c.items()})
    
        # update total likelihood and f vector
        n_reads = len(L_r_c)
        log_L_r = {r:math.log(v)-log_c[r] for r,v in L_r_c.items()}    
        prev_log_likelihood = total_log_likelihood
        total_log_likelihood = sum(log_L_r.values())
        f = {s:(sum(r_map.values())/n_reads) for s,r_map in L_sgr.items()}
        
        print("*****ITERATION*****")
        print(total_log_likelihood)
        
        # check f vector sums to 1
        f_sum = sum(f.values())
        if not (.999 <= f_sum <= 1.0001):
            raise ValueError(f"f sums to {f_sum}, rather than 1")
    
        # confirm log likelihood increase
        log_likelihood_diff = total_log_likelihood - prev_log_likelihood
        if log_likelihood_diff < 0:
            raise ValueError(f"total_log_likelihood decreased from prior iteration")
    
        # exit loop if small increase
        if total_log_likelihood - prev_log_likelihood < .001:
            return f
        
def f_reduce(f, threshold):
    '''reduce composition vector f to only those with value > threshold, then normalize
    
        f: composition vector 
        threshold: float to cut off value
        return: array of reduced composition vector
    '''
    f_dropped = {k:v for k, v in f.items() if v > threshold}
    f_total = sum(f_dropped .values())
    return {k:v/f_total for k,v in f_dropped.items()}

def f_to_lineage_df(f, seq2taxid_path, tsv_output_name):
    '''converts composition vector f to a pandas df with cols lineages, abundance.
        Stores df as .tsv file in tsv_output_name
    
        f (arr): composition vector with keys database sequnce names for abundace output
        seq2taxid_path (str): string path to tab separated file of seqid and taxids
        tsv_output_name (str): name of output .tsv file of generated dataframe
        returns (df): pandas df with lineage and abundances for values in f
    '''
    ncbi = NCBITaxa()
    df_seq2tax = pd.read_csv(seq2taxid_path, sep='\t', header=None)
    seq2taxid = dict(zip(df_seq2tax[0], df_seq2tax[1])) 
    ncbi.translate_to_names(seq2taxid.values())

    lineages = [(';').join(ncbi.translate_to_names(ncbi.get_lineage(seq2taxid[k]))) for k in f]
    results_df = pd.DataFrame(list(zip(lineages, f.values())), 
               columns =['lineage', 'abundance'])
    results_df = results_df[results_df['abundance']>0].groupby("lineage").sum().sort_values('lineage')
    results_df.to_csv(f"{tsv_output_name}.tsv", sep='\t')
    return results_df


def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
                    'input_file', type=str,
                    help='filepath to input [fasta,fastq,sam]')
    args = parser.parse_args()


    # initialize values
    db_fasta_path = "ncbi16s_db/bacteria_and_archaea.16SrRNA.fna"
    db_seq2taxid_path = "ncbi16s_db/ncbi16s_seq2tax.map"
    db_ids = [record.id for record in SeqIO.parse(db_fasta_path, "fasta")]
    output_dir = "results/"
    if not os.path.exists(output_dir): os.makedirs(output_dir)
    
    filename = pathlib.PurePath(args.input_file).stem 
    filetype = pathlib.PurePath(args.input_file).suffix
    pwd = os.getcwd()
    
    if filetype == '.sam':
        sam_file = f"{args.input_file}"
    else:
        sam_file = os.path.join(output_dir, f"{filename}.sam")
        pwd = os.getcwd()
        subprocess.check_output(
            f"bwa mem -a -x ont2d {pwd}/ncbi16s_db/bwa_dbs/ncbi_16s {args.input_file} > {sam_file}",
            shell=True)


    ## script
    p_char = get_char_align_probabilites(sam_file)
    log_L_rgs = log_L_rgs_dict(sam_file, p_char)
    f = EM_iterations(log_L_rgs, db_ids)
    f_dropped = f_reduce(f, .01)
    results_df = f_to_lineage_df(f_dropped, db_seq2taxid_path, os.path.join(output_dir, filename))
    
main()



