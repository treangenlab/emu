## EMU: relative abundance estimator for full-length 16s reads


### Description

Emu is a realtive abundance estimator..

### Synopsis

- Calculate abundances Oxford Nanopore single-end:
```bash
python EM_classification.py example/full_length.fa
```
- Calculate abundances short-read:
```bash
python EM_classification.py --sr example/short_read1.fa example/short_read2.fa
```

### Installation

bioconda...

```bash
conda config --add channels r
conda config --add channels bioconda
conda install pysam
conda install pandas
pip install flatten-dict

```

### Options





### Build Custom Database

An emu database consists of 4 files:

- names_df.tsv: tab separated datasheet of database sequence names containing at least columns: 'tax_id' and 'name_txt'
- nodes_df.tsv: tab separated datasheet of database sequence lineages containing at least columns: 'tax_id', 'parent_tax_id', and 'rank'
- species_taxid.fasta: database sequences where each sequence header starts with the repsective species-level tax_id preceeding a colon [\<species_taxid>:\<remainder of header>]
- unique_taxids.tsv: single column tab separated values of unqiue tax_ids in database

To build a custom database with corresponding NCBI taxonomy, 4 files are needed.

- names.dmp: names file from NCBI taxonomy dump
- nodes.dmp: nodes file from NCBI taxonomy dump
- database.fasta: nucleotide sequences
- seq2taxid.map: headerless two column tab-separated values, where each row contains (1) sequence header in database.fasta and (2) species-level tax id.

```bash
python build_database.py <db_name> --names <names.dmp> --nodes <nodes.dmp> --sequences <database.fasta> --seq2tax <seq2taxid.map>
```

Example:

```bash
python build_database.py zymo_assembled_db --names example_customdb/names.dmp --nodes example_customdb/nodes.dmp --sequences ./example_customdb/zymo_assembled.fasta --seq2tax ./example_customdb/zymo_assembled_seq2tax.map
```

```bash
python EM_classification.py example/emu_example.fa --db ./zymo_assembled_db
```


### Limitations




