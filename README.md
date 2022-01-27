## Emu: species-level taxonomic abundance for full-length 16S reads


### Description

Emu is a relative abundance estimator for 16S genomic sequences. The method is optimized for error-prone full-length reads, but can also be utilized for short-read data.

### Demo

Calculate relative abundances for Oxford Nanopore Technologies single-end 16S reads:
```bash
emu abundance example/full_length.fa
```
Calculate relative abundances for short paired-end 16S data:
```bash
emu abundance --type sr example/short_read_f.fq example/short_read_r.fq
```
Calculate relative abundances for short single-end 16S data:
```bash
emu abundance --type sr example/short_read_f.fq
```

Expected output: each of the above commands is expected to create a relative abundance .tsv file in a ./results folder. Output names are: full_length_rel-abundance.tsv, short_read_f-short_read_r_rel-abundance.tsv, short_read_f_rel-abundance.tsv.

Expected run time: each of the above commands is expected to complete successfully in a couple seconds.
### Installation

##### 1. Download database

Define `<path_to_database>` as your desired database directory. If you desire a different database, skip this step and follow steps below in **Build Custom Database**.

```bash
export EMU_DATABASE_DIR=<path_to_database>
wget -qO- https://gitlab.com/treangenlab/emu/-/archive/v1.0.1/emu-v3.0.0.tar.gz | tar -C $EMU_DATABASE_DIR -xvz --strip-components=2 emu-v3.0.0/emu_database/
```

##### 2. Activate appropriate conda environment

Emu requires Python version to be in range 3.6.0-3.8.12. 

###### Option A: Create new Conda environment

```bash
conda create --name py37 python=3.7 
conda activate py37
```

###### Option B: Set Python version in current environment

```bash
conda install python=3.7
```

##### 3. Install Emu

###### Option A: Install Emu via conda

Add the bioconda channel and install Emu.
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install emu
```

###### Option B: Create local Emu conda environment

If you are unable to install Emu via conda as described above, an alternative approach is to install conda and create an environment that supports Emu. The default name of the environment created is `emu`, but this can be configured in the `environment.yml` file if desired. The environment will need to be activated before Emu can be run.

```bash
conda env create -f environment.yml
conda activate custom-emu
```
NOTE: with this installation method, all commands will need to be run with `<path to emu>` instead of `emu`. For example, in the directory containing the Emu script, a working command would be:
```bash
./emu abundance example/full_length.fa
```

Each step of the installation process is expected to take a matter of seconds.

### Abundance Estimation Parameters

| Command	| Default	| Description	|
| :-------  | :----- | :-------- | 
|--type	| map-ont	| denote sequencer [short-read:sr, Pac-Bio:map-pb, ONT:map-ont]	|
|--min-abundance| 0.0001| generates results with species relative abundance above this value in addition to full results; .01 = 1%|
|--db| $EMU_DATABASE_DIR| path to emu database; directory must include the following files: names_df.tsv, nodes_df.tsv, species_taxid.fasta, unqiue_taxids.tsv|
|--N| 50| max number of alignments utilized for each read|
|--output-dir| ./results| directory for output results|
|--output-basename| stem of input_file(s)| basename of all output files saved in output-dir; default utilizes basename from input file(s)|
|--keep-files| FALSE| keep working files in output-dir ( alignments [.sam], reads of specied length [.fa])|
|--threads| 3| number of threads utilized by minimap2|

Note: If you are experiencing heavy RAM consumption, first upgrade minimap2 to at least v2.22. If memory is still an issue, try decreasing the number of secondary alignments evaluated for each read (--N).

### Build Custom Database

An emu database consists of 4 files:
| Filename	| Description	|
| :-------  | :----- |
|names_df.tsv| tab separated datasheet of database taxonomy names containing at least columns: 'tax_id' and 'name_txt'|
|nodes_df.tsv| tab separated datasheet of database taxonomy tree containing at least columns: 'tax_id', 'parent_tax_id', and 'rank'|
|species_taxid.fasta| database sequences where each sequence header starts with the respective species-level tax id (or lowest level above species-level if missing) preceeding a colon [&lt;species_taxid>:&lt;remainder of header>]|
|unique_taxids.tsv| single column tab separated values of unqiue tax_ids in database|

To build a custom database with corresponding NCBI taxonomy, 4 files are needed.

- names.dmp: names file from NCBI taxonomy dump
- nodes.dmp: nodes file from NCBI taxonomy dump
- database.fasta: nucleotide sequences
- seq2taxid.map: headerless two column tab-separated values, where each row contains (1) sequence header in database.fasta and (2) tax id.

```bash
emu build-database <db_name> --names <names.dmp> --nodes <nodes.dmp> --sequences <database.fasta> --seq2tax <seq2taxid.map>
```

Example:

```bash
emu build-database zymo_assembled_db --names ./example_customdb/ex_names.dmp --nodes ./example_customdb/ex_nodes.dmp --sequences ./example_customdb/ex.fasta --seq2tax ./example_customdb/ex_seq2tax.map
```

```bash
emu abundance ./example_customdb/ex.fasta --db ./zymo_assembled_db
```

If preferred, user can define custom database through the shell variable rather than specifying with each run at the command line.

```bash
export EMU_DATABASE_DIR=./zymo_assembled_db
emu abundance ./example_customdb/ex.fasta
```

Note: If your taxonomy is missing species-level information, a “pseudo” species will be reported as “unclassified &lt;genus>” where &lt;genus> is the labeled genus in the taxonomic lineage. If genus-level classification is also missing in the lineage, this process will continue moving up the taxonomic lineage until a specified label (&lt;taxa>) is detected. Then, "unclassified &lt;taxa>" will be reported as the species classification instead. 

### Collapse Taxonomy

The collapse-taxonomy function can be used on any emu output &lt;.tsv> file to generate an additional output collapsed at the desired taxonomic rank. File is output the same folder as the input file, with filename:&lt;input_file>-&lt;rank>.tsv. Accepted ranks: ['species', 'genus', 'family', 'order', 'class', 'phylum', 'clade', 'superkingdom']

```bash
emu collapse-taxonomy <file_path> <rank>
```

### System Requirements

All software depencies are listed in environment.yml. Emu v2.0.1 has been tested on Python v3.7 and used to generate results in manuscript.

### Emu Manuscript

Preprint: [Curry, Kristen, et al. "Emu: Species-Level Microbial Community Profiling for Full-Length Nanopore 16S Reads." bioRxiv (2021).](https://www.biorxiv.org/content/10.1101/2021.05.02.442339v1)
Repository for reproduction of results in manuscript: [Emu-benchmark](https://gitlab.com/treangenlab/emu-benchmark)



