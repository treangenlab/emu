#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kcurry
"""

import os
import argparse
import pathlib
import subprocess
from sys import stdout
from operator import add, mul

import math
import pysam
import numpy as np
import pandas as pd
from flatten_dict import unflatten

## static global variables
CIGAR_OPS = [1, 4, 8]   #['I', 'S', 'X']


def get_align_stats(alignment):
    """Retrieve list of inquired cigar stats for alignment ['I', 'S', '=', 'X']

        alignment (pysam.AlignmentFile): align of interest (requires =/X CIGAR operators)
        return (list(int)): list of counts for each cigar operation defined in CIGAR_OPS
    """
    return [alignment.get_cigar_stats()[0][cigar_op] for cigar_op in CIGAR_OPS]


def get_cigar_op_log_probabilites(sam_path):
    """P(match), P(mismatch), P(insertion), P(softclipping),
        by counting how often the corresponding operations occur in the primary alignments
        and by normalizing over the total sum of operations

        sam_path(str): path to sam file of interest (requires =/X CIGAR operators)
        return: log probabilites (list(float)) for each cigar operation defined in CIGAR_OPS,
                    where p > 0
                zero_locs (list(int)): list of indices (int) where probability == 0
    """
    cigar_stats_primary = [0] * len(CIGAR_OPS)
    # pylint: disable=maybe-no-member
    sam_pysam = pysam.AlignmentFile(sam_path)
    for alignment in sam_pysam.fetch():
        if not alignment.is_secondary and not alignment.is_supplementary \
                and alignment.reference_name:
            cigar_stats_primary = list(map(add, cigar_stats_primary, get_align_stats(alignment)))
    zero_locs = [i for i, e in enumerate(cigar_stats_primary) if e == 0]

    if zero_locs:
        for i in sorted(zero_locs, reverse=True):
            del cigar_stats_primary[i]
    n_char = sum(cigar_stats_primary)
    return [math.log(x) for x in np.array(cigar_stats_primary)/n_char], zero_locs


def compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op):
    """ log(L(r|s)) = log(P(cigar_op)) × n_cigar_op for cigar_ops: ['I', 'S', '=', 'X']

        alignment(pysam.AlignmentFile): pysam alignment to score
        cigar_stats(list(int)): list of cigar stats to compute
        log_p_cigar_op(list(float)): list of cigar_op probabilities corresponding to cigar_stats;
                                        computed from primary alignments
        return: log_score (float): log(L(r|s))
                query_name (str): query name in alignment
                species_tid (int): species-level taxonomy id corresponding to ref_name
    """

    log_score = sum(list(map(mul, log_p_cigar_op, cigar_stats)))
    ref_name, query_name = alignment.reference_name, alignment.query_name
    species_tid = int(ref_name.split(":")[0])
    return log_score, query_name, species_tid


def log_prob_rgs_dict(sam_path, log_p_cigar_op, p_cigar_op_zero_locs=None):
    """dict containing log(L(read|seq)) for all pairwise alignments in sam file

        sam_path(str): path to sam file
        log_p_cigar_op(list(float)): probability for each cigar operation defined in CIGAR_OPS,
                                        where p > 0
        zero_locs(list(int)): list of indices (int) where probability == 0
        return ({[str,int]:float}): dict[(query_name,ref_tax_id)]=log(L(query_name|ref_tax_id))
    """
    # calculate log(L(read|seq)) for all alignments
    log_p_rgs = {}
    # pylint: disable=maybe-no-member
    sam_filename = pysam.AlignmentFile(sam_path, 'rb')
    unassigned_count = 0

    if not p_cigar_op_zero_locs:
        for alignment in sam_filename.fetch():
            if alignment.reference_name:
                cigar_stats = get_align_stats(alignment)
                log_score, query_name, species_tid = \
                    compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op)
                if (((query_name, species_tid) not in log_p_rgs or
                     log_p_rgs[(query_name, species_tid)] < log_score)):
                    log_p_rgs[(query_name, species_tid)] = log_score
            else:
                unassigned_count += 1
    else:
        for alignment in sam_filename.fetch():
            if alignment.reference_name:
                cigar_stats = get_align_stats(alignment)
                if sum([cigar_stats[x] for x in p_cigar_op_zero_locs]) == 0:
                    for i in sorted(p_cigar_op_zero_locs, reverse=True):
                        del cigar_stats[i]
                    log_score, query_name, species_tid = \
                        compute_log_prob_rgs(alignment, cigar_stats, log_p_cigar_op)
                    if (((query_name, species_tid) not in log_p_rgs or
                         log_p_rgs[(query_name, species_tid)] < log_score)):
                        log_p_rgs[(query_name, species_tid)] = log_score
            else:
                unassigned_count += 1

    stdout.write(f"Unassigned read count: {unassigned_count}\n")
    return log_p_rgs


def expectation_maximization(log_p_rgs, freq):
    """One iteration of the EM algorithm. Updates the relative abundance estimation in f based on
            probabolities in log_L_rgs.

        log_p_rgs({[str,int]:float}): dict[(query_name,ref_tax_id)]=log(L(query_name|ref_tax_id))
        freq{int:float}: dict[species_tax_id]:likelihood species is present in sample
        returns: f {int:float}: dict[species_tax_id]:updated likelihood species is present in sample
                total_log_likelihood (float): log likelihood updated f is accurate
    """
    # log_p_rns = log(L(r|s))+log(f(s)) for each alignment
    log_p_rns = {}
    for (read, seq) in log_p_rgs:
        if seq in freq and freq[seq] != 0:
            log_p_rns[(read, seq)] = (log_p_rgs[(read, seq)] + math.log(freq[seq]))
    log_c = {r: -max(smap.values()) for r, smap in unflatten(log_p_rns).items()}
    p_rns_c = {(r, s): (math.e ** (v + log_c[r])) for (r, s), v in log_p_rns.items()}
    p_r_c = {r: sum(s.values()) for r, s in unflatten(p_rns_c).items()}

    # P(s|r), likelihood read r emanates db seq s
    p_sgr = unflatten({(seq, read): v / p_r_c[read] for (read, seq), v in p_rns_c.items()})

    # update total likelihood and f vector
    log_p_r = {read: math.log(v) - log_c[read] for read, v in p_r_c.items()}
    n_reads = len(log_p_r)
    return {seq: (sum(r_map.values()) / n_reads) for seq, r_map in p_sgr.items()},\
           sum(log_p_r.values())


def expectation_maximization_iterations(log_p_rgs, db_ids, lli_thresh, input_threshold):
    """Full expectation maximization algorithm for alignments in log_L_rgs dict

        log_p_rgs{[str,int]:float}: dict[(query_name,ref_tax_id)]=log(L(query_name|ref_tax_id))
        db_ids(list(int)): list of each unique species taxonomy id present in database
        lli_thresh(float): log likelihood increase minimum to continue EM iterations
        input_threshold(float): minimum relative abundance in output
        return: {int:float}: dict[species_tax_id]:estimated likelihood species is present in sample
    """
    n_db = len(db_ids)
    n_reads = len(unflatten(log_p_rgs))
    stdout.write(f"Assigned read count: {n_reads}\n")
    freq, counter = dict.fromkeys(db_ids, 1 / n_db), 1

    # set output abundance threshold
    freq_thresh = 1/n_reads
    if n_reads > 1000:
        freq_thresh = 10/n_reads

    total_log_likelihood = -math.inf
    while True:
        freq, updated_log_likelihood = expectation_maximization(log_p_rgs, freq)

        # check f vector sums to 1
        freq_sum = sum(freq.values())
        if not .9 <= freq_sum <= 1.1:
            raise ValueError(f"f sums to {freq_sum}, rather than 1")

        # confirm log likelihood increase
        log_likelihood_diff = updated_log_likelihood - total_log_likelihood
        total_log_likelihood = updated_log_likelihood
        if log_likelihood_diff < 0:
            raise ValueError("total_log_likelihood decreased from prior iteration")

        # exit loop if log likelihood increase less than threshold
        if log_likelihood_diff < lli_thresh:
            stdout.write(f"Number of EM iterations: {counter}\n")
            freq = {k: v for k, v in freq.items() if v >= freq_thresh}
            freq_full, updated_log_likelihood = expectation_maximization(log_p_rgs, freq)
            freq_set_thresh = None
            if freq_thresh < input_threshold:
                freq = {k: v for k, v in freq_full.items() if v >= input_threshold}
                freq_set_thresh, updated_log_likelihood = expectation_maximization(log_p_rgs, freq)
            return freq_full, freq_set_thresh

        counter += 1

def lineage_dict_from_tid(taxid, nodes_df, names_df):
    """Retrieve dict of lineage for given taxid

        tid(int): tax id to retrieve lineage dict
        nodes_df(df): pandas df of nodes.dmp with columns ['tax_id', 'parent_tax_id', 'rank'];
                            tax_id as index
        names_df(df): pandas df of names.dmp with columns ['tax_id', 'name_txt']; tax_id as index
        return {str:str}: dict[rank]:taxonomy name at rank
    """

    lineage_dict = {}
    rank = None
    while rank != "no rank":
        row = nodes_df.loc[taxid]
        rank = row["rank"]
        lineage_dict[rank] = names_df.loc[taxid]["name_txt"]
        taxid = row["parent_tax_id"]
    return lineage_dict


def freq_to_lineage_df(freq, tsv_output_path, nodes_df, names_df):
    """Converts freq to a pandas df where each row contains abundance and tax lineage for
                classified species in f.keys(). Stores df as .tsv file in tsv_output_path.

        freq{int:float}: dict[species_tax_id]:estimated likelihood species is present in sample
        tsv_output_path(str): path and name of output .tsv file for generated dataframe
        nodes_df(df): pandas df of nodes.dmp with columns ['tax_id', 'parent_tax_id', 'rank'];
                            tax_id as index
        names_df(df): pandas df of names.dmp with columns ['tax_id', 'name_txt']; tax_id as index
        returns(df): pandas df with lineage and abundances for values in f
    """

    results_df = pd.DataFrame(list(zip(freq.keys(), freq.values())),
                              columns=["tax_id", "abundance"])
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
        'input_file', type=str, nargs='+',
        help='filepath to input nt sequence file')
    parser.add_argument(
        '--short-read', '-s', action='store_true',
        help='apply tag if short read data')
    parser.add_argument(
        '--min-read-len', type=int, default=1000,
        help='minimum read length for long-read data')
    parser.add_argument(
        '--max-read-len', type=int, default=5000,
        help='maximum read length for long-read data')
    parser.add_argument(
        '--min-abundance', '-a', type=float, default=0.0001,
        help='min species abundance in results [0.0001]')
    parser.add_argument(
        '--db', type=str, default="./emu_database",
        help='path to emu database containing: names_df.tsv, '
             'nodes_df.tsv, species_taxid.fasta, unqiue_taxids.tsv')
    parser.add_argument(
        '--N', type=int, default=25,
        help='minimap max number of alignments per read [25]')
    parser.add_argument(
        '--output-dir', type=str, default="./",
        help='output directory name [./]')
    parser.add_argument(
        '--output', '-o', type=str,
        help='output filename [{input_file}]')
    parser.add_argument(
        '--threads', type=int, default=3,
        help='threads utilized by minimap [3]')
    args = parser.parse_args()

    # convert taxonomy files to dataframes
    df_nodes = pd.read_csv(os.path.join(args.db, "nodes_df.tsv"), sep='\t').set_index('tax_id')
    df_names = pd.read_csv(os.path.join(args.db, "names_df.tsv"), sep='\t').set_index('tax_id')
    db_species_tids = pd.read_csv(os.path.join(args.db, "unique_taxids.tsv"), sep='\t')['tax_id']

    # set up output paths
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    out_file = os.path.join(args.output_dir, f"{pathlib.PurePath(args.input_file[0]).stem}")
    if args.output:
        out_file = os.path.join(args.output_dir, args.output)
    filetype = pathlib.PurePath(args.input_file[0]).suffix

    # generate read alignments
    INPUT_FILE = " ".join(args.input_file)
    if filetype == '.sam':
        sam_file = f"{INPUT_FILE}"
    else:
        sam_file = f"{out_file}.sam"
        db_sequence_file = os.path.join(args.db, 'species_taxid.fasta')
        if args.short_read:
            subprocess.check_output(
                f"minimap2 -x sr -ac -t {args.threads} -N {args.N} "
                f"-p .9 --eqx {db_sequence_file} {INPUT_FILE} -o {sam_file}",
                shell=True)
        else:
            fasta_trimmed = f"{out_file}_trimmed.fa"
            subprocess.check_output(
                f"bioawk -c fastx '(length($seq)<={args.max_read_len}"
                f"&&length($seq)>={args.min_read_len}){{print \">\" $name ORS $seq}}' "
                f"{INPUT_FILE} > {fasta_trimmed}",
                shell=True)
            subprocess.check_output(
                f"minimap2 -x map-ont -ac -t {args.threads} -N {args.N} -p .9 --eqx "
                f"{db_sequence_file} {fasta_trimmed} -o {sam_file}",
                shell=True)
            os.remove(fasta_trimmed)

    # perform EM algorithm & generate output
    log_prob_cigar_op, locs_p_cigar_zero = get_cigar_op_log_probabilites(sam_file)
    log_prob_rgs = log_prob_rgs_dict(sam_file, log_prob_cigar_op, locs_p_cigar_zero)
    f_full, f_set_thresh = expectation_maximization_iterations\
        (log_prob_rgs, db_species_tids, .01, args.min_abundance)
    results_df_full = freq_to_lineage_df(f_full, out_file, df_nodes, df_names)
    if f_set_thresh:
        results_df_thresh = freq_to_lineage_df(f_set_thresh,
                                            f"{out_file}_abundance_thresh_{args.min_abundance}",
                                            df_nodes, df_names)

    # clean up extra file
    if os.path.exists(sam_file):
        os.remove(sam_file)
