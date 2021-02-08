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
import numpy as np
import pandas as pd
from flatten_dict import unflatten
from operator import add, mul
from Bio import SeqIO

## static global variables
CIGAR_OPS = [1, 4, 7, 8]


def get_align_stats(alignment):
    """Convert list of CIGAR stats containing NM to a dict containing mismatches
        
        alignment: pysam.AlignmentFile
        return: list of counts (ints) for each cigar operation defined in CIGAR_OPS
    """
    return [alignment.get_cigar_stats()[0][cigar_op] for cigar_op in CIGAR_OPS]


def get_cigar_op_log_probabilites(sam):
    """P(match), P(mismatch), P(insertion), P(deletion), P(clipping),
        by counting how often the corresponding operations occur in the primary alignments
        and by normalizing over the total sum of operations
    
        sam: path to sam file bwa output
        return: list of log probabilites (float) for each cigar operation defined in CIGAR_OPS, where p > 0
                zero_locs: list of indices (int) where probability == 0
    """
    cigar_stats_primary = [0] * len(CIGAR_OPS)
    samfile = pysam.AlignmentFile(sam)
    for alignment in samfile.fetch():
        if not alignment.is_secondary and not alignment.is_supplementary and alignment.reference_name:
            cigar_stats_primary = list(map(add, cigar_stats_primary, get_align_stats(alignment)))
    zero_locs = [i for i, e in enumerate(cigar_stats_primary) if e == 0]

    if zero_locs:
        for i in sorted(zero_locs, reverse=True):
            del cigar_stats_primary[i]
    n_char = sum(cigar_stats_primary)
    return [math.log(x) for x in np.array(cigar_stats_primary)/n_char], zero_locs


def compute_log_L_rgs(alignment, cigar_stats, log_p_cigar_op):
    """ log(L(r|s)) = log(P(match)) × n_matches + log(P(mismatches)) × n_mismatches ...
    
        alignment: pysam alignment to score
        cigar_stats: list of cigar stats to compute
        log_p_cigar_op: list of probabilities values corresponding to cigar_stats
        return: log_score (float), ref_name (str), query_name (str), specied_tid (int)
    """

    log_score = sum(list(map(mul, log_p_cigar_op, cigar_stats)))
    ref_name, query_name = alignment.reference_name, alignment.query_name
    species_tid = int(ref_name.split(":")[0])
    return log_score, ref_name, query_name, species_tid


def log_L_rgs_dict(sam_path, log_p_cigar_op, p_cigar_op_zero_locs=None):
    """dict containing log(L(r|s)) for all pairwise alignments in bwa output
    
        sam_path: path to sam file
        p_char: dict of likelihood for each cigar output with prob > 0
        return: dict[(r,s)]=log(L(r|s))
    """
    # calculate log(L(r|s)) for all alignments
    log_L_rgs = {}
    samfile = pysam.AlignmentFile(sam_path)

    if not p_cigar_op_zero_locs:
        for alignment in samfile.fetch():
            if alignment.reference_name:
                cigar_stats = get_align_stats(alignment)
                log_score, ref_name, query_name, species_tid = compute_log_L_rgs(alignment, cigar_stats, log_p_cigar_op)
                if (((query_name, species_tid) not in log_L_rgs or log_L_rgs[(query_name, species_tid)] < log_score)):
                    log_L_rgs[(query_name, species_tid)] = log_score
    else:
        for alignment in samfile.fetch():
            if alignment.reference_name:
                cigar_stats = get_align_stats(alignment)
                if sum([cigar_stats[x] for x in p_cigar_op_zero_locs]) == 0:
                    for i in sorted(p_cigar_op_zero_locs, reverse=True):
                        del cigar_stats[i]
                    log_score, ref_name, query_name, species_tid = compute_log_L_rgs(alignment, cigar_stats, log_p_cigar_op)
                    if (((query_name, species_tid) not in log_L_rgs or log_L_rgs[(query_name, species_tid)] < log_score)):
                        log_L_rgs[(query_name, species_tid)] = log_score
    return log_L_rgs


def EM_iterations(log_L_rgs, db_ids, lli_thresh, names_df, nodes_df, input_threshold, output_dir, fname):
    """Expectation maximization algorithm for alignments in log_L_rgs dict
        
        log_L_rgs: dict[(r,s)]=log(L(r|s))
        db_ids: array of ids in database
        lli_thresh: log likelihood increase minimum to continue EM iterations
        names_df: pandas df of names.dmp with columns ['tax_id', 'name_txt','unique_name', 'name_class'] with tax_id as index
        nodes_df: pandas df of nodes.dmp with columns ['tax_id', 'parent_tax_id', 'rank'] with tax_id as index
        input_threshold: float of minimum relative abundance in output
        output_dir: str path of output directory
        fname: str output filename

        return: composition vector f[s]=relative abundance of s from sequence database
    """
    n_db = len(db_ids)
    n_reads = len(unflatten(log_L_rgs))
    f = dict.fromkeys(db_ids, 1 / n_db)
    counter, break_flag = 1, False

    # set up dir to output results after each iteration
    dir = f"{os.path.join(output_dir, fname)}_iterations"
    if not os.path.exists(dir):
        os.makedirs(dir)

    # set output abundance threshold
    f_thresh = 1/n_reads
    if n_reads > 1000:
        f_thresh = 10/n_reads
    if input_threshold < f_thresh:
        f_thresh = input_threshold

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
        log_L_r = {r: math.log(v) - log_c[r] for r, v in L_r_c.items()}
        prev_log_likelihood = total_log_likelihood
        total_log_likelihood = sum(log_L_r.values())
        f = {s: (sum(r_map.values()) / n_reads) for s, r_map in L_sgr.items()}

        print(f"*****ITERATION:{counter}*****")
        print(total_log_likelihood)

        # check f vector sums to 1
        f_sum = sum(f.values())
        if not (.999 <= f_sum <= 1.0001):
            raise ValueError(f"f sums to {f_sum}, rather than 1")

        if break_flag:
            print(f"Number of EM iterations: {counter}")
            return f

        # confirm log likelihood increase
        log_likelihood_diff = total_log_likelihood - prev_log_likelihood
        if log_likelihood_diff < 0:
            raise ValueError(f"total_log_likelihood decreased from prior iteration")

        # exit loop if log likelihood increase less than threshold
        if log_likelihood_diff < lli_thresh:
            f = {k: v for k, v in f.items() if v >= f_thresh}
            break_flag = True

        f_to_lineage_df(f, f"{dir}/{counter}", nodes_df, names_df)
        counter += 1


def create_nodes_df(nodes_path):
    """convert nodes.dmp file into pandas dataframe

        nodes_path (str): path to nodes.dmp file
        returns (df): pandas df of nodes.dmp with columns ['tax_id', 'parent_tax_id', 'rank'] with tax_id as index
    """
    node_headers = ['tax_id', 'parent_tax_id', 'rank']
    nodes_df = pd.read_csv(nodes_path, sep='\t', header=None, dtype=str)[[0, 2, 4]]
    nodes_df.columns = node_headers
    nodes_df = nodes_df.set_index("tax_id")
    return nodes_df


def create_names_df(names_path):
    """convert names.dmp file into pandas dataframe

        names_path (str): path to names.dmp file
        returns (df): pandas df of names.dmp with columns ['tax_id', 'name_txt','unique_name', 'name_class'] with tax_id as index
    """
    name_headers = ['tax_id', 'name_txt', 'unique_name', 'name_class']
    names_df = pd.read_csv(names_path, sep='\t', index_col=False, header=None, dtype=str).drop([1, 3, 5, 7], axis=1)
    names_df.columns = name_headers
    names_df = names_df[names_df["name_class"] == "scientific name"].set_index("tax_id")
    return names_df


def get_fasta_ids(fasta_path):
    """ Returns list of unique ids in fasta [fasta_path]

        fasta_path (str): path to database fasta
        returns (list): list of unique ids (str) in fasta
    """
    tids = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        tids += [record.id.split(":")[0]]
    return np.unique(tids)


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


def f_to_lineage_df(f, tsv_output_path, nodes_df, names_df):
    """converts composition vector f to a pandas df where each row contains abundance and
        tax lineage for each classified species.
        Stores df as .tsv file in tsv_output_name.

        f (arr): composition vector with keys database sequnce names for abundace output
        tsv_output_path (str): path and name of output .tsv file of generated dataframe
        nodes_df (df): pandas df of nodes.dmp with columns ['tax_id', 'parent_tax_id', 'rank'] with tax_id as index
        names_df (df): pandas df of names.dmp with columns ['tax_id', 'name_txt','unique_name', 'name_class'] with tax_id as index
        returns (df): pandas df with lineage and abundances for values in f
    """

    results_df = pd.DataFrame(list(zip(f.keys(), f.values())), columns=["tax_id", "abundance"])
    lineages = results_df["tax_id"].apply(lambda x: lineage_dict_from_tid(x, nodes_df, names_df))
    results_df = pd.concat([results_df, pd.json_normalize(lineages)], axis=1)
    header_order = ["abundance", "species", "genus", "family", "order", "class",
                    "phylum", "clade", "superkingdom", "subspecies",
                    "species subgroup", "species group", "tax_id"]
    for col in header_order:
        if col not in results_df.columns:
            results_df[col] = ""
    results_df = results_df.sort_values(header_order[8:0:-1]).reset_index(drop=True)
    results_df = results_df.reindex(header_order, axis=1)

    results_df.to_csv(f"{tsv_output_path}.tsv", index=False, sep='\t')
    return results_df


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input_file', type=str,
        help='filepath to input [fasta,fastq,sam]')
    parser.add_argument(
        '--lli', '-l', type=float, default=0.1,
        help='min log likelihood increase to continue EM iterations [0.1]')
    parser.add_argument(
        '--threshold', '-t', type=float, default=0.0001,
        help='min species abundance in results [0.0001]')
    parser.add_argument(
        '--names', type=str, default="./database/NCBI_taxonomy/names.dmp",
        help='path to names.dmp')
    parser.add_argument(
        '--nodes', type=str, default="./database/NCBI_taxonomy/nodes.dmp",
        help='path to nodes.dmp')
    parser.add_argument(
        '--threads', type=int, default=40,
        help='threads utilized by minimap')
    parser.add_argument(
        '--db', type=str, default="./database/combined_tid.fasta",
        help='path to fasta file of database sequences')
    parser.add_argument(
        '--output', '-o', type=str,
        help='output filename')
    parser.add_argument(
        '--output_dir', type=str, default="results_test/",
        help='output directory name')
    args = parser.parse_args()

    # convert taxonomy files to dataframes
    database = 'default'
    if database == 'default':
        nodes_df = pd.read_csv("./database/NCBI_taxonomy/nodes_df.tsv", sep='\t').set_index('tax_id')
        names_df = pd.read_csv("./database/NCBI_taxonomy/names_df.tsv", sep='\t').set_index('tax_id')
        db_species_tids = pd.read_csv("./database/unique_taxids.tsv", sep='\t')['taxonomy_id']
    else: ##change how input database in handled
        nodes_df = create_nodes_df(args.nodes)
        names_df = create_names_df(args.names)
        db_species_tids = get_fasta_ids(args.db)


    # output files
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    filename = f"{pathlib.PurePath(args.input_file).stem}"
    if args.output:
        filename = args.output
    filetype = pathlib.PurePath(args.input_file).suffix
    pwd = os.getcwd()

    # input files
    if filetype == '.sam':
        sam_file = f"{args.input_file}"
    else:
        sam_file = os.path.join(args.output_dir, f"{filename}.sam")
        pwd = os.getcwd()
        subprocess.check_output(
            f"minimap2 -x map-ont -ac -t {args.threads} -N 50 -p .9 --eqx {args.db} {args.input_file} -o {sam_file}",
            shell=True)

    # script
    log_p_cigar_op, p_cigar_zero_locs = get_cigar_op_log_probabilites(sam_file)
    log_L_rgs = log_L_rgs_dict(sam_file, log_p_cigar_op, p_cigar_zero_locs)
    f = EM_iterations(log_L_rgs, db_species_tids, args.lli, names_df, nodes_df, args.threshold, args.output_dir, filename)
    results_df_full = f_to_lineage_df(f, f"{os.path.join(args.output_dir, filename)}", nodes_df, names_df)
    #results_df_full = f_to_lineage_df_lineage_txt(f, f"{os.path.join(args.output_dir, filename)}_full_{args.lli}")
    #results_df = df_reduce(results_df_full, args.threshold)
    #results_df.to_csv(f"{os.path.join(args.output_dir, filename)}.tsv", index=False, sep='\t')
