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
import re

import math
import pysam
import pandas as pd
from collections import Counter
from flatten_dict import unflatten
from Bio import SeqIO


def start_clip_length(align):
    """align: pysam alignment
    
        return: int of starting clip length
    """
    if align.cigartuples[0][0] in [4, 5]:
        return align.cigartuples[0][1]
    return 0


def get_n_locs(align):
    """align: pysam alignment
        
        return: list of positions in query which align to ambiguous base
    """
    n_cols = []
    if align.has_tag("cs"):
        i = start_clip_length(align)
        cs_list = re.findall(r'([+*:-]{1})(\w+)', align.get_tag("cs"))
        for (k, v) in cs_list:
            if k == ":":
                i += int(v)
            elif k == "+":
                i += len(v)
            elif k == "*":
                if "n" in v: n_cols += [i]
                i += 1
    return n_cols


def create_n_cols_query_dict(sam):
    """sam: path to samfile
    
        return: dict[query] = set of positions which align to ambiguous base
    """
    samfile = pysam.AlignmentFile(sam)
    n_cols_dict = {}
    for align in samfile.fetch():
        if align.query_name not in n_cols_dict:
            n_cols_dict[align.query_name] = set()
        if align.reference_name:
            n_match_list = get_n_locs(align)
            if n_match_list:
                n_cols_dict[align.query_name] = n_cols_dict[align.query_name].union(set(n_match_list))
    return n_cols_dict


def cigar_to_remove(align, n_cols_dict):
    """align: pysam alignment
        n_cols_dict: dict[query] = set of positions which align to ambiguous base
        
        return: dict of cigar stats for alignment at positions in n_cols_dict
    """
    cols_remove = sorted(n_cols_dict[align.query_name]) + [None]
    col_to_remove = cols_remove.pop(0)
    cigar_removed = dict.fromkeys([":", "+", "*"], 0)
    if align.has_tag("cs"):
        i = start_clip_length(align)
        cs_list = re.findall(r'([+*:-]{1})(\w+)', align.get_tag("cs"))
        for (k, v) in cs_list:
            if not col_to_remove:
                break
            if k == ":":
                i += int(v)
            elif k == "+":
                i += len(v)
            elif k == "*":
                i += 1
            while col_to_remove < i:
                cigar_removed[k] += 1
                col_to_remove = cols_remove.pop(0)
                if not col_to_remove: break
    cs_to_cigar_map = {"+": "I", "*": "X", ":": "="}
    for k, v in cs_to_cigar_map.items():
        cigar_removed[v] = cigar_removed.pop(k)
    return cigar_removed


def get_align_stats(alignment):
    """Convert list of CIGAR stats containing NM to a dict containing mismatches
        
        alignment: pysam.AlignmentFile
        return: dict of cigar stats
    """
    cigar_stats = alignment.get_cigar_stats()[0]
    cigar_headers = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B', 'NM']
    char_align_stats = dict(zip(cigar_headers, cigar_stats))
    char_align_stats['C'] = char_align_stats['H'] + char_align_stats['S']
    for k in ['M', 'N', 'S', 'H', 'P', 'B', 'NM']:
        char_align_stats.pop(k)
    return char_align_stats


def get_char_align_probabilites(sam, remove_n_cols=False):
    """P(match), P(mismatch), P(insertion), P(deletion), P(clipping),
        by counting how often the corresponding operations occur in the primary alignments
        and by normalizing over the total sum of operations
    
        sam: path to sam file bwa output
        return: dict of likelihood for each cigar output with prob > 0
    """
    # initialize character alignment probabilities
    remove_cols_dict = None
    if remove_n_cols:
        remove_cols_dict = create_n_cols_query_dict(sam)
    readlist, char_align_stats_primary = [], {}
    samfile = pysam.AlignmentFile(sam)
    for read in samfile.fetch():
        if not read.is_secondary and not read.is_supplementary:
            if read.query_name in readlist:
                print("ERROR!!")
            readlist += [read.query_name]
            read_align_stats = get_align_stats(read)
            if remove_n_cols:
                cigar_removed_stats = cigar_to_remove(read, remove_cols_dict)
                for k, v in cigar_removed_stats.items():
                    read_align_stats[k] -= v
            char_align_stats_primary = Counter(char_align_stats_primary) + Counter(read_align_stats)
    char_align_stats_primary_oi = {k: char_align_stats_primary[k] for k in ['=', 'C', 'X', 'I']}
    n_char = sum(char_align_stats_primary_oi.values())
    return {char: val / n_char for char, val in char_align_stats_primary_oi.items()}, remove_cols_dict


def compute_log_L_rgs(p_char, cigar_stats):
    """ log(L(r|s)) = log(P(match)) ×n_matches + log(P(mismatches)) ×n_mismatches ...
    
        p_dict: dict of CIGAR likelihoods
        cigar_stats: dict of cigar stats to compute
        return: float log(L(r|s)) for sequences r and s in cigar_stats alignment
    """
    cigar_stats_oi = {k: cigar_stats[k] for k in ['=', 'C', 'X', 'I']}
    #log_char_sum = math.log(sum(cigar_stats_oi.values()))
    char_sum = sum(cigar_stats_oi.values())
    value = 0
    prob_sum = 0
    for char, count in cigar_stats_oi.items():
        if char not in p_char and cigar_stats_oi[char] > 0:
            return None
        elif char in p_char:
            prob_sum += p_char[char] * cigar_stats_oi[char]
            value += math.log(p_char[char]) * cigar_stats_oi[char]
    #value += math.log(prob_sum/sum(cigar_stats_oi.values())) * cigar_stats['I']
    return value


def log_L_rgs_dict(bwa_sam, p_char, remove_cols_dict=None):
    """dict containing log(L(r|s)) for all pairwise alignments in bwa output
    
        bwa_sam: path to sam file bwa output
        p_char: dict of likelihood for each cigar output with prob > 0
        return: dict[(r,s)]=log(L(r|s))
    """
    # calculate log(L(r|s)) for all alignments
    log_L_rgs = {}
    samfile = pysam.AlignmentFile(bwa_sam)
    data_cigars = []
    for alignment in samfile.fetch():
        if alignment.reference_name:
            align_stats = get_align_stats(alignment)
            data = [alignment.query_name, alignment.reference_name, align_stats.copy()]
            if remove_cols_dict:
                cigar_removed_stats = cigar_to_remove(alignment, remove_cols_dict)
                for k, v in cigar_removed_stats.items():
                    align_stats[k] = align_stats[k] - v
            val = compute_log_L_rgs(p_char, align_stats)
            data += [align_stats, val]
            data_cigars.append(data)
            if ((alignment.query_name, alignment.reference_name) not in log_L_rgs
                or log_L_rgs[(alignment.query_name, alignment.reference_name)] < val) and val:
                log_L_rgs[(alignment.query_name, alignment.reference_name)] = val
        else:
            data_cigars.append([alignment.query_name, "-"])

    df_cigars = pd.DataFrame(data_cigars, columns=['query', 'reference', 'cigar', 'cigar_n_removed', 'score'])
    df_cigars.to_csv("cigar_scores.tsv", sep='\t', index=False)
    return log_L_rgs


def EM_iterations(log_L_rgs, db_ids):
    """Expectation maximization algorithm for alignments in log_L_rgs dict
        
        log_L_rgs: dict[(r,s)]=log(L(r|s))
        db_ids: array of ids in database
        n_reads: int number of reads in input sequence (|R|)
        return: composition vector f[s]=relative abundance of s from sequence database 
    """
    n_db = len(db_ids)
    f = dict.fromkeys(db_ids, 1 / n_db)

    total_log_likelihood = -math.inf
    while (True):
        # log_L_rns = log(L(r|s))+log(f(s)) for each alignment
        log_L_rns = {}
        for (r, s) in log_L_rgs:
            if s in f and f[s] != 0:
                log_L_rns[(r, s)] = (log_L_rgs[(r, s)] + math.log(f[s]))
        log_c = {r: -max(smap.values()) for r, smap in unflatten(log_L_rns).items()}
        L_rns_c = {(r, s): (math.e ** (v + log_c[r])) for (r, s), v in log_L_rns.items()}
        L_r_c = {r: sum(s.values()) for r, s in unflatten(L_rns_c).items()}

        # P(s|r), likelihood read r emanates db seq s
        L_sgr = unflatten({(s, r): v / L_r_c[r] for (r, s), v in L_rns_c.items()})

        # update total likelihood and f vector
        n_reads = len(L_r_c)
        log_L_r = {r: math.log(v) - log_c[r] for r, v in L_r_c.items()}
        prev_log_likelihood = total_log_likelihood
        total_log_likelihood = sum(log_L_r.values())
        f = {s: (sum(r_map.values()) / n_reads) for s, r_map in L_sgr.items()}

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
        if total_log_likelihood - prev_log_likelihood < 1:
            return f


def f_reduce(f, threshold):
    """reduce composition vector f to only those with value > threshold, then normalize
    
        f: composition vector 
        threshold: float to cut off value
        return: array of reduced composition vector
    """
    f_dropped = {k: v for k, v in f.items() if v > threshold}
    f_total = sum(f_dropped.values())
    return {k: v / f_total for k, v in f_dropped.items()}


def lineage_dict_from_tid(tid, nodes_df, names_df):
    """get dict of lineage for given taxid
    
    tid(str): tax id to extract lineage dict from
    nodes_df(df): pandas df of nodes.dmp with columns ['tax_id', 'parent_tax_id', 'rank'] with tax_id as index
    names_df(df): pandas df of names.dmp with columns ['tax_id', 'name_txt','unique_name', 'name_class'] with tax_id as index
    return(dict): [taxonomy rank]:name
    """
    lineage_dict = {}
    rank = None
    while rank != "superkingdom":
        row = nodes_df.loc[tid]
        rank = row["rank"]
        lineage_dict[rank] = names_df.loc[tid]["name_txt"]
        tid = row["parent_tax_id"]
    return lineage_dict


def f_to_lineage_df(f, tsv_output_name, nodes_path, names_path, seq2taxid_path):
    """converts composition vector f to a pandas df where each row contains abundance and
        tax lineage for each classified species.
        Stores df as .tsv file in tsv_output_name.

        f (arr): composition vector with keys database sequnce names for abundace output
        tsv_output_name (str): name of output .tsv file of generated dataframe
        nodes_path (str): string path to nodes.dmp from ncbi taxonomy dump
        names_path (str): string path to names.dmp from ncbi taxonomy dump
        seq2taxid_path (str): string path to tab separated file of seqid and taxids
        returns (df): pandas df with lineage and abundances for values in f
    """
    name_headers = ['tax_id', 'name_txt', 'unique_name', 'name_class']
    node_headers = ['tax_id', 'parent_tax_id', 'rank']

    # convert taxonomy files to dataframes
    seq2tax_df = pd.read_csv(seq2taxid_path, sep='\t', header=None)
    seq2taxid = dict(zip(seq2tax_df[0], seq2tax_df[1]))
    names_df = pd.read_csv(names_path, sep='\t', index_col=False, header=None, dtype=str).drop([1, 3, 5, 7], axis=1)
    names_df.columns = name_headers
    names_df = names_df[names_df["name_class"] == "scientific name"].set_index("tax_id")
    nodes_df = pd.read_csv(nodes_path, sep='\t', header=None, dtype=str)[[0, 2, 4]]
    nodes_df.columns = node_headers
    nodes_df = nodes_df.set_index("tax_id")

    tax_ids = [seq2taxid[seqid] for seqid in f.keys()]
    f_tax = dict.fromkeys(tax_ids, 0)
    for seqid, v in f.items():
        f_tax[seq2taxid[seqid]] += v
    results_df = pd.DataFrame(list(zip(f_tax.keys(), f_tax.values())), columns=["tax_id", "abundance"])
    lineages = results_df["tax_id"].apply(lambda x: lineage_dict_from_tid(str(x), nodes_df, names_df))
    results_df = pd.concat([results_df, pd.json_normalize(lineages)], axis=1)
    header_order = ["abundance", "species", "genus", "family", "order", "class",
                    "phylum", "clade", "superkingdom", "strain", "subspecies",
                    "species subgroup", "species group", "tax_id"]
    results_df = results_df.sort_values(header_order[8:0:-1]).reset_index(drop=True)
    results_df = results_df.reindex(header_order, axis=1)

    results_df.to_csv(f"{tsv_output_name}.tsv", index=False, sep='\t')
    return results_df


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input_file', type=str,
        help='filepath to input [fasta,fastq,sam]')
    args = parser.parse_args()

    # initialize values
    db_tax_path = "ncbi16s_db/"
    names_path = db_tax_path + "NCBI_taxonomy/names.dmp"
    nodes_path = db_tax_path + "NCBI_taxonomy/nodes.dmp"
    #seq2taxid_path = db_tax_path + "ncbi16s_seq2tax.map"
    seq2taxid_path = "input/ncbi16s_zymo_db/ncbi16s_seq2tax.map"
    #db_fasta_path = "ncbi16s_db/bacteria_and_archaea.16SrRNA.fna"
    db_fasta_path = "input/ncbi16s_zymo_db/bac_arch_zymo.16SrRNA.fna"
    db_ids = [record.id for record in SeqIO.parse(db_fasta_path, "fasta")]
    output_dir = "results_new/"
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
            f"minimap2 -x map-ont -ac -t 40 -N 100 --eqx --cs ncbi16s_db/bacteria_and_archaea.16SrRNA.fna {args.input_file} -o {sam_file}",
            # f"bwa mem -a -x ont2d {pwd}/ncbi16s_db/bwa_dbs/ncbi_16s {args.input_file} > {sam_file}",
            shell=True)

    # script
    p_char, remove_cols_dict = get_char_align_probabilites(sam_file, remove_n_cols=False)
    log_L_rgs = log_L_rgs_dict(sam_file, p_char, remove_cols_dict)
    f = EM_iterations(log_L_rgs, db_ids)
    f_dropped = f_reduce(f, .01)
    results_df_full = f_to_lineage_df(f, f"{os.path.join(output_dir, filename)}_full", nodes_path, names_path,
                                      seq2taxid_path)
    results_df = f_to_lineage_df(f_dropped, os.path.join(output_dir, filename), nodes_path, names_path, seq2taxid_path)
