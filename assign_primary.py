
import os
import subprocess
import argparse
import pathlib
import pandas as pd

import pysam
from Bio import SeqIO
import EM_classification as utils


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
    seq2taxid_path = "input/ncbi16s_zymo_db/ncbi16s_seq2tax.map"
    seq2tax_df = pd.read_csv(seq2taxid_path, sep='\t', header=None)
    seq2taxid = dict(zip(seq2tax_df[0], seq2tax_df[1]))
    seq2taxid['-'] = '-'
    names_df = pd.read_csv(names_path, sep='\t', index_col=False, header=None, dtype=str).drop([1, 3, 5, 7], axis=1)
    names_df = names_df[names_df[6] == "scientific name"]
    tax2name = dict(zip(names_df[0], names_df[2]))
    tax2name['-'] = '-'
    #db_fasta_path = "input/ncbi16s_zymo_db/zymo_assmebled_only.16SrRNA.fna"
    db_fasta_path = "ncbi16s_db/bacteria_and_archaea.16SrRNA.fna"
    db_ids = [record.id for record in SeqIO.parse(db_fasta_path, "fasta")]
    output_dir = "results/"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    filename = pathlib.PurePath(args.input_file).stem
    filetype = pathlib.PurePath(args.input_file).suffix
    pwd = os.getcwd()

    if filetype == '.sam':
        sam_file = f"{args.input_file}"
    else:
        sam_file = os.path.join(output_dir, f"{filename}.sam")
        pwd = os.getcwd()
        subprocess.check_output(
            f"minimap2 -x map-ont -ac -t 40 -N 100 --eqx --cs {db_fasta_path} {args.input_file} -o {sam_file}",
            # f"bwa mem -a -x ont2d {pwd}/ncbi16s_db/bwa_dbs/ncbi_16s {args.input_file} > {sam_file}",
            shell=True)

    # initialize
    n_db = len(db_ids)
    f_counts = dict.fromkeys(db_ids, 0)
    assignments = []

    # add one count for each primary alignment
    samfile = pysam.AlignmentFile(sam_file)
    for align in samfile.fetch():
        if not align.is_secondary and not align.is_supplementary and align.reference_name:
            f_counts[align.reference_name] = f_counts[align.reference_name] + 1
            assignments.append([align.query_name, align.reference_name])
        if not align.reference_name:
            assignments.append([align.query_name, "-"])

    # convert counts to percentage
    total = sum(f_counts.values())
    f = {k: (v/total) for k, v in f_counts.items() if v > 0}
    df_assignments = pd.DataFrame(assignments, columns=["query", "reference"])
    df_assignments["tax_id"] = df_assignments["reference"].apply(lambda x: seq2taxid[x])
    df_assignments["species_name"] = df_assignments["tax_id"].apply(lambda x: tax2name[str(x)])


    results_df_full = utils.f_to_lineage_df(f, f"{os.path.join(output_dir, filename)}_primary_full", nodes_path, names_path,
                                      seq2taxid_path)
    df_assignments.to_csv(f"{os.path.join(output_dir, filename)}_primary_assignments.tsv", sep='\t', index=False)
