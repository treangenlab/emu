#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kcurry
"""

import os
from sys import stdout

import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def create_nodes_df(nodes_path):
    """convert nodes.dmp file into pandas dataframe

        nodes_path(str): path to nodes.dmp file
        returns(df): pandas df containing 'tax_id', 'parent_tax_id', and 'rank'
    """
    node_headers = ['tax_id', 'parent_tax_id', 'rank']
    nodes_df = pd.read_csv(nodes_path, sep='\t', header=None, dtype=str)[[0, 2, 4]]
    nodes_df.columns = node_headers
    nodes_df = nodes_df.set_index("tax_id")
    return nodes_df


def create_names_df(names_path):
    """convert names.dmp file into pandas dataframe

        names_path(str): path to names.dmp file
        returns(df): pandas df containing ['tax_id','name_txt'] for each row with
                    'name_class' = 'scientific name'
    """
    name_headers = ['tax_id', 'name_txt', 'unique_name', 'name_class']
    names_df = pd.read_csv(names_path, sep='\t', index_col=False, header=None, dtype=str)\
        .drop([1, 3, 5, 7], axis=1)
    names_df.columns = name_headers
    names_df = names_df[names_df["name_class"] == "scientific name"].set_index("tax_id")
    names_df = names_df.drop(columns=['unique_name', 'name_class'])
    return names_df

def get_species_tid(tid, nodes_df):
    """ Get species level taxid in lineage for taxid [tid]

        tid(int): taxid for species level or more specific
        nodes_df(pandas df): pandas df containing columns ['tax_id', 'parent_tax_id', 'rank']
        return(int): species taxid in lineage
    """
    if str(tid) not in nodes_df.index:
        raise ValueError(f"Taxid:{tid} not found in nodes file")
    row = nodes_df.loc[str(tid)]
    while row['rank'] != 'species':
        parent_tid = row['parent_tax_id']
        row = nodes_df.loc[str(parent_tid)]
        if row['rank'] == 'superkingdom':
            return None
    return row.name

def create_seq2tax_dict(seq2tax_path, nodes_df):
    """Convert seqid-taxid mapping in seq2tax_path to dict mapping seqid to species level taxid

        seq2tax_path(str): path to seqid-taxid mapping file
        nodes_df(df): pandas df containing columns ['tax_id', 'parent_tax_id', 'rank']
        returns {str:int}: dict[seqid] = species taxid
    """
    seq2tax_dict, species_id_dict = {}, {}
    with open(seq2tax_path) as file:
        for line in file:
            (seqid, tid) = line.rstrip().split("\t")
            if tid in species_id_dict.keys():
                species_tid = species_id_dict[tid]
            else:
                species_tid = get_species_tid(tid, nodes_df)
                species_id_dict[tid] = species_tid
            seq2tax_dict[seqid] = species_tid
    return seq2tax_dict

def create_unique_seq_dict(db_fasta_path, seq2tax_dict):
    """ Creates dict of unique sequences to species taxids connected with the sequence

        db_fasta_path(str): path to fasta file of database sequences
        seq2tax_dict{str:int}: dict[seqid] = species taxid
        returns {str:{int:[str]}}: dict[seq] = {species_taxid: [list of sequence ids]}
    """
    fasta_dict = {}
    for record in SeqIO.parse(db_fasta_path, "fasta"):
        tid = seq2tax_dict[record.id]
        if tid:
            if record.seq in fasta_dict.keys():
                if tid in fasta_dict[record.seq].keys():
                    fasta_dict[record.seq][tid] += [record.description]
                else:
                    fasta_dict[record.seq][tid] = [record.description]
            elif record.seq.reverse_complement() in fasta_dict.keys():
                if tid in fasta_dict[record.seq.reverse_complement()].keys():
                    fasta_dict[record.seq.reverse_complement()][tid] += [record.description]
                else:
                    fasta_dict[record.seq.reverse_complement()][tid] = [record.description]
            else:
                fasta_dict[record.seq] = {tid: [record.description]}
    return fasta_dict

def create_reduced_fasta(fasta_dict, db_name):
    """ Creates fasta file of taxid for each sequences in fasta_dict with id
            'species_taxid:db_name:sequence_id'

        fasta_dict{str:{int:[str]}}: dict[seq] = {species_taxid: [list of sequence ids]}
        db_name(str): name to represent database represented in fasta_dict
        returns (list[Bio.SeqRecord]): list of sequences for output fasta file
    """
    records, count = [], 1
    for seq, tid_dict in fasta_dict.items():
        for taxid, descriptions in tid_dict.items():
            records += [SeqRecord(seq,
                                  id=f"{taxid}:{db_name}:{count}", description=f"{descriptions}")]
            count += 1
    return records

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'db_name', type=str,
        help='custom database name')
    parser.add_argument(
        '--names', type=str,
        help='path to names.dmp file')
    parser.add_argument(
        '--nodes', type=str,
        help='path to nodes.dmp file')
    parser.add_argument(
        '--sequences', type=str,
        help='path to fasta of database sequences')
    parser.add_argument(
        '--seq2tax', type=str,
        help='path to tsv mapping species tax id to fasta sequence headers')
    args = parser.parse_args()

    emu_db_path = os.getcwd()
    custom_db_path = os.path.join(emu_db_path, args.db_name)
    if not os.path.exists(custom_db_path):
        os.makedirs(custom_db_path)
    stdout.write(f"Emu custom database generating at path: {custom_db_path} ...\n")

    df_names = create_names_df(args.names)
    df_nodes = create_nodes_df(args.nodes)
    seq2tax = create_seq2tax_dict(args.seq2tax, df_nodes)
    dict_fasta = create_unique_seq_dict(args.sequences, seq2tax)
    db_unique_ids = set(seq2tax.values())
    fasta_records = create_reduced_fasta(dict_fasta, args.db_name)

    df_names.to_csv(os.path.join(custom_db_path, 'names_df.tsv'), sep='\t')
    df_nodes.to_csv(os.path.join(custom_db_path, 'nodes_df.tsv'), sep='\t')
    pd.DataFrame(db_unique_ids, columns=['tax_id']).\
        to_csv(os.path.join(custom_db_path, 'unique_taxids.tsv'),index=False, sep='\t')
    SeqIO.write(fasta_records, os.path.join(custom_db_path, 'species_taxid.fasta'), "fasta")
    stdout.write("Database creation successful\n")
