from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd

def get_species_tid(tid, nodes_df):
    if str(tid) not in nodes_df.index:
        print(tid)
    else:
        row = nodes_df.loc[str(tid)]
        while row['rank'] != 'species':
            parent_tid = row['parent_tax_id']
            row = nodes_df.loc[str(parent_tid)]
        return row.name

def seq2tax_dict(seq2tax_path, nodes_df):
    seq2tax = {}
    with open(seq2tax_path) as f:
        for line in f:
            (seqid, tid) = line.rstrip().split("\t")
            species_tid = get_species_tid(tid, nodes_df)
            seq2tax[seqid] = species_tid
    return seq2tax


def create_nodes_df(nodes_path):
    node_headers = ['tax_id', 'parent_tax_id', 'rank']
    nodes_df = pd.read_csv(nodes_path, sep='\t', header=None, dtype=str)[[0, 2, 4]]
    nodes_df.columns = node_headers
    nodes_df = nodes_df.set_index("tax_id")
    return nodes_df

def create_unique_seq_dict(db_fasta_path, seq2tax_dict):
    fasta_dict = {}
    for record in SeqIO.parse(db_fasta_path, "fasta"):
        #name = seq2tax_dict[record.description.split("|")[1]]
        name = seq2tax_dict[record.id]
        if record.seq in fasta_dict.keys():
            if name in fasta_dict[record.seq].keys():
                fasta_dict[record.seq][name] += 1
            else:
                fasta_dict[record.seq][name] = 1
        elif record.seq.reverse_complement() in fasta_dict.keys():
            if name in fasta_dict[record.seq.reverse_complement()].keys():
                fasta_dict[record.seq.reverse_complement()][name] += 1
            else:
                fasta_dict[record.seq.reverse_complement()][name] = 1
        else:
            fasta_dict[record.seq] = {name: 1}
    return fasta_dict

def create_reduced_fasta(fasta_dict, db_name):
    records = []
    rec_num = 0
    for k, d in fasta_dict.items():
        if len(d) == 1:
            rec_num += 1
            taxid = list(d.keys())[0]
            records += [SeqRecord(k, id=f"{taxid}:{db_name}:{rec_num}", description=taxid)]
    SeqIO.write(records, f"{db_name}_speciesid_removedups.fasta", "fasta")


db_name = "curated"
db_fasta_path = "db_curated/curated_db.fasta"
#db_seq_to_taxid_map = "rrndb/seq2taxid.map"
#nodes_path = "ncbi16s_db/NCBI_taxonomy/nodes.dmp"

#nodes_df = create_nodes_df(nodes_path)
#seq2tax_dict = seq2tax_dict(db_seq_to_taxid_map, nodes_df)
fasta_dict = create_unique_seq_dict(db_fasta_path, seq2tax_dict=None)
create_reduced_fasta(fasta_dict, db_name)


diff_spec_dict = {}
for k, d in fasta_dict.items():
    if len(d) > 1:
        diff_spec_dict[k] = d

with open(f'sheet_multi_species_seqs_{db_name}.csv', 'w') as f:
    for k in diff_spec_dict.keys():
        f.write("%s,%s\n" % (k, diff_spec_dict[k]))



