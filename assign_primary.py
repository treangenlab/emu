
import os
import subprocess
import argparse
import pathlib
import pandas as pd

import pysam
import EM_classification as utils


def f_to_lineage_df(f, tsv_output_name, nodes_df, names_df):
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

    results_df = pd.DataFrame(list(zip(f.keys(), f.values())), columns=["tax_id", "counts"])
    lineages = results_df["tax_id"].apply(lambda x: utils.lineage_dict_from_tid(str(x), nodes_df, names_df))
    results_df = pd.concat([results_df, pd.json_normalize(lineages)], axis=1)
    header_order = ["counts", "species", "genus", "family", "order", "class",
                    "phylum", "clade", "superkingdom", "strain", "subspecies",
                    "species subgroup", "species group", "tax_id"]
    for col in header_order:
        if col not in results_df.columns:
            results_df[col] = ""
    results_df = results_df.sort_values(header_order[8:0:-1]).reset_index(drop=True)
    results_df = results_df.reindex(header_order, axis=1)

    return results_df


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        'input_file', type=str,
        help='filepath to input [sam]')
    args = parser.parse_args()

    # initialize values
    names_path = "./NCBI_taxonomy/names.dmp"
    nodes_path = "./NCBI_taxonomy/nodes.dmp"
    names_df = pd.read_csv(names_path, sep='\t', index_col=False, header=None, dtype=str).drop([1, 3, 5, 7], axis=1)
    names_df = names_df[names_df[6] == "scientific name"]
    tax2name = dict(zip(names_df[0], names_df[2]))
    tax2name['-'] = '-'
    output_dir = "./"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # convert taxonomy files to dataframes
    name_headers = ['tax_id', 'name_txt', 'unique_name', 'name_class']
    node_headers = ['tax_id', 'parent_tax_id', 'rank']
    names_df = pd.read_csv(names_path, sep='\t', index_col=False, header=None, dtype=str).drop([1, 3, 5, 7], axis=1)
    names_df.columns = name_headers
    names_df = names_df[names_df["name_class"] == "scientific name"].set_index("tax_id")
    nodes_df = pd.read_csv(nodes_path, sep='\t', header=None, dtype=str)[[0, 2, 4]]
    nodes_df.columns = node_headers
    nodes_df = nodes_df.set_index("tax_id")

    filename = pathlib.PurePath(args.input_file).stem
    filetype = pathlib.PurePath(args.input_file).suffix
    pwd = os.getcwd()

    if filetype != '.sam':
        raise ValueError("Input file must be of type <.sam>")

    sam_file = f"{args.input_file}"

    # initialize
    f_counts, assignments = {}, []
    unassigned_count = 0

    # add one count for each primary alignment
    samfile = pysam.AlignmentFile(sam_file)
    for align in samfile.fetch():
        if not align.reference_name:
            assignments.append([align.query_name, "-"])
            unassigned_count += 1
        elif not align.is_secondary and not align.is_supplementary:
            tid = align.reference_name.split(":")[0]
            if tid not in f_counts.keys():
                f_counts[tid] = 1
            else:
                f_counts[tid] = f_counts[tid] + 1
            assignments.append([align.query_name, align.reference_name])

    # convert counts to percentage
    total = sum(f_counts.values())
    total_with_unassigned = total + unassigned_count
    #f = {k: (v/total) for k, v in f_counts.items() if v > 0}
    results_df_full = f_to_lineage_df(f_counts, f"{os.path.join(output_dir, filename)}_primary_full", nodes_df, names_df)
    results_df_full['abundance'] = results_df_full['counts'].apply(lambda x: x/total)
    results_df_full = results_df_full.append(pd.DataFrame({"counts": [unassigned_count],
                                                           "species": ["unclassified"],
                                                           "genus":["unclassified"]}), ignore_index=True)
    results_df_full['abundance with unassigned'] = results_df_full['counts'].apply(lambda x: x / total_with_unassigned)
    results_df_full.to_csv(f"{os.path.join(output_dir, filename)}_truth.tsv", sep='\t', index=False)

    df_assignments = pd.DataFrame(assignments, columns=["query", "reference"])
    df_assignments["tax_id"] = df_assignments["reference"].apply(lambda x: x.split(":")[0])
    df_assignments["species_name"] = df_assignments["tax_id"].apply(lambda x: tax2name[str(x)])
    df_assignments.to_csv(f"{os.path.join(output_dir, filename)}_primary_assignments.tsv", sep='\t', index=False)
