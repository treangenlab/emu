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
wget -qO- https://gitlab.com/treangenlab/emu/-/archive/v3.0.0/emu-v3.0.0.tar.gz | tar -C $EMU_DATABASE_DIR -xvz --strip-components=2 emu-v3.0.0/emu_database/
```

** Note Emu v3.0+ database requirements differ from previous versions. Check you are using the appropriate database for the version you are running. Both databases contain identical information: a combination of [rrnDB v5.6](https://rrndb-umms-med-umich-edu.translate.goog/?_x_tr_sl=en&_x_tr_tl=fr&_x_tr_hl=fr&_x_tr_pto=sc) and [NCBI 16S RefSeq](https://www.ncbi.nlm.nih.gov/refseq/targetedloci/16S_process/) from 17 September, 2020. Taxonomy is also from NCBI on the same date. The resulting database contains 49,301 sequences from 17,555 unique bacterial and archaeal species.

##### 2. Activate appropriate conda environment

Emu requires Python version to be >=3.6. 

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
Although not required, we recommend using minimap2 version >=2.22 to avoid memory leak [bug with highly repetitive sequences](https://github.com/lh3/minimap2/releases).

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
|--N| 50| max number of alignments utilized for each read in minimap2|
|--K| 500M| minibatch size for mapping in minimap2|
|--output-dir| ./results| directory for output results|
|--output-basename| stem of input_file(s)| basename of all output files saved in output-dir; default utilizes basename from input file(s)|
|--keep-files| FALSE| keep working files in output-dir ( alignments [.sam], reads of specied length [.fa])|
|--keep-counts| FALSE| include estimated read counts for each species in output|
|--keep-read-assignments| FALSE| output .tsv file with read assignment distributions: each row as an input read; each entry as the likelihood it is dervied from that taxa (taxid is the column header); each row sums to 1|
|--output-unclassified| FALSE| generate a separate output file of unclassified sequences|
|--threads| 3| number of threads utilized by minimap2|


Note: If you are experiencing heavy RAM consumption, first upgrade minimap2 to at least v2.22. If memory is still an issue, try decreasing the number of secondary alignments evaluated for each read (--N).
Note: Estimated read counts are based on likelihood probabilities and therefore may not be integer values.

### Build Custom Database

An emu database consists of 2 files:
| Filename	| Description	|
| :-------  | :----- |
|taxonomy.tsv| tab separated datasheet of database taxonomy lineages containing at columns: 'tax_id' and any taxonomic ranks (i.e. species, genus, etc) |
|species_taxid.fasta| database sequences where each sequence header starts with the respective species-level tax id (or lowest level above species-level if missing) preceeding a colon [&lt;species_taxid>:&lt;remainder of header>]|

A custom database can be built from either NCBI taxonomy (names.dmp and nodes.dmp files) or a .tsv file containing full taxonomy lineages.
If the database if build from NCBI taxonomy, the database fasta sequences will be classified at the species level (or the lowest taxonomic rank
provided if already classified above the species level). If taxonomy is provided as a list of lineage, database fasta sequences
will be classified at the rank of the tax_id provided in the seq2taxid.map file. Therefore, when using direct taxonomy to build
your custom database, we recommend providing a seq2taxid.map at the species-level or higher.

The following files are required to build a custom database:
| Command	| file(s)	| Description	|
| :-------  | :----- | :-------- | 
|--sequences	| database.fasta	| nucleotide sequences	|
|--seq2tax	| database.fasta	| headerless two column tab-separated values, where each row contains (1) sequence header in database.fasta and (2) corresponding tax id	|
| *either taxonomy option:*| 
|--ncbi-taxonomy	| names.dmp & nodes.dmp	| directory containing both names.dmp and nodes.dmp files in NCBI taxonomy format and named accordingly	|
|--taxonomy-list	| input taxonomy.tsv	| a .tsv file containing complete taxonomic lineages. The first column MUST be the taxonomy ids. Remaining columns can be in any format, then Emu abundance output will match this format|


```bash
emu build-database <db_name> --sequences <database.fasta> --seq2tax <seq2taxid.map> --ncbi-taxonomy <dir-to-names/nodes.dmp>
OR
emu build-database <db_name> --sequences <database.fasta> --seq2tax <seq2taxid.map> --taxonomy-list <taxonomy.tsv>
```

Example:

```bash
emu build-database zymo_assembled_db --sequences ./example_customdb/ex.fasta --seq2tax ./example_customdb/ex_seq2tax.map --ncbi-taxonomy ./example_customdb/
OR
emu build-database zymo_assembled_db --sequences ./example_customdb/ex.fasta --seq2tax ./example_customdb/ex_seq2tax.map --taxonomy-list ./example_customdb/ex-taxonomy.tsv
```

```bash
emu abundance ./example_customdb/ex.fasta --db ./zymo_assembled_db --threads 3
```

If preferred, user can define custom database through the shell variable rather than specifying with each run at the command line.

```bash
export EMU_DATABASE_DIR=./zymo_assembled_db
emu abundance ./example_customdb/ex.fasta
```

Note for NCBI-taxonomy created database: If your taxonomy is missing species-level information, a “pseudo” species will be reported as “unclassified &lt;genus>” where &lt;genus> is the labeled genus in the taxonomic lineage. If genus-level classification is also missing in the lineage, this process will continue moving up the taxonomic lineage until a specified label (&lt;taxa>) is detected. Then, "unclassified &lt;taxa>" will be reported as the species classification instead. 

#### Alternative Databases

[RDP v11.5](https://rdp.cme.msu.edu/) with NCBI taxonomy has been pre-built for Emu v3.0+ and can be downloaded accordingly. 

```bash
export EMU_DATABASE_DIR=<path_to_database>
wget -qO- https://gitlab.com/treangenlab/emu/-/archive/v3.0.0/emu-v3.0.0.tar.gz | tar -C $EMU_DATABASE_DIR -xvz --strip-components=2 emu-v3.0.0/rdp_database/
```


[SILVA v138.1](https://www.arb-silva.de/) has been pre-built for Emu v3.0+ from the [DADA2 SILVA species-level database](https://zenodo.org/record/4587955#.YvqmSezMLOQ).
```bash
export EMU_DATABASE_DIR=<path_to_database>
wget -qO- https://gitlab.com/treangenlab/emu/-/archive/v3.4.1/emu-v3.4.1.tar.gz | tar -C $EMU_DATABASE_DIR -xvz --strip-components=2 emu-v3.4.1/silva_database/
```

[UNITE general fasta v8.3 fungi](https://plutof.ut.ee/#/doi/10.15156/BIO/1280049) has been pre-built for Emu v3.0+. 
This database has not yet been tested or validated with Emu.
```bash
export EMU_DATABASE_DIR=<path_to_database>
wget -qO- https://gitlab.com/treangenlab/emu/-/archive/v3.4.2/emu-v3.4.2.tar.gz | tar -C $EMU_DATABASE_DIR -xvz --strip-components=2 emu-v3.4.2/unite-fungi_database/
```

[UNITE general fasta v8.3 all eukaryotes](https://plutof.ut.ee/#/doi/10.15156/BIO/1280127) has been pre-built for Emu v3.0+. 
This database has not yet been tested or validated with Emu.
```bash
export EMU_DATABASE_DIR=<path_to_database>
wget -qO- https://gitlab.com/treangenlab/emu/-/archive/v3.4.2/emu-v3.4.2.tar.gz | tar -C $EMU_DATABASE_DIR -xvz --strip-components=2 emu-v3.4.2/unite-all_database/
```

### Collapse Taxonomy

The collapse-taxonomy function can be used on any emu output &lt;.tsv> file to generate an additional output collapsed 
at the desired taxonomic rank. File is output the same folder as the input file, with filename:&lt;input_file>-&lt;rank>.tsv. 
Accepted ranks: ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']

```bash
emu collapse-taxonomy <file_path> <rank>
```

### Combine Outputs

The combine-outputs function can be used to create a single table containing all Emu output relative abundances in a single directory.
Note this function will select all the .tsv files in the provided directory that contain 'rel-abundance' in the filename.
Combined table will only include all ranks above the specified rank according to this list:
[tax_id, species, genus, family, order, class, phylum, superkingdom]. Specified rank must
be in this list and in each of the included Emu outputs. Combined table will be created in the provided directory path with the file name:
emu-combined-&lt;rank>.tsv. In order to include tax_id in your output, specific &lt;rank> as "tax_id".

```bash
emu combine-outputs <directory_path> <rank>
```

Optional additional parameters:
| Command	| Description	|
| :-------  | :-------- | 
|--split-tables	| output 2 tables: (1) abundances only at specified rank and (2) taxonomic lineages down to specified rank	|
|--counts	| output estimated counts rather than relative abundance percentage in combined table. Only includes Emu relative abundance outputs that already have 'estimated counts'  |

### System Requirements

All software dependencies are listed in environment.yml. Emu v3.0.0 has been tested on Python v3.7 and used to generate results in manuscript.

### Emu Manuscript

Publication: [Kristen D. Curry et al., “Emu: Species-Level Microbial Community Profiling of Full-Length 16S RRNA Oxford Nanopore Sequencing Data,” Nature Methods, June 30, 2022, 1–9, https://doi.org/10.1038/s41592-022-01520-4.](https://www.nature.com/articles/s41592-022-01520-4)
Repository for reproduction of results in manuscript: [Emu-benchmark](https://gitlab.com/treangenlab/emu-benchmark)

### Database Citations
Please use citations below if any of the pre-contructed databases are utilized:

##### Emu default database
- Stoddard, S. F., Smith, B. J., Hein, R., Roller, B. R. & Schmidt, T. M. (2015) rrnDB: improved tools for interpreting rRNA gene abundance in bacteria and archaea and a new foundation for future development. Nucl. Acids Res. 43, D593–D598.
- O’Leary, N. A. et al. (2016) Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucl. Acids Res. 44, D733–D745.
- Schoch, C. L. et al. (2020) NCBI Taxonomy: a comprehensive update on curation, resources and tools. Database (Oxford) 2020.


##### RDP
- Cole, J. R. et al. (2014) Ribosomal Database Project: data and tools for high throughput rRNA analysis. Nucl. Acids Res. 42, D633–D642.
- Schoch, C. L. et al. (2020) NCBI Taxonomy: a comprehensive update on curation, resources and tools. Database (Oxford) 2020.

##### SILVA
- Quast C, Pruesse E, Yilmaz P, Gerken J, Schweer T, Yarza P, Peplies J, Glöckner FO (2013) The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. Nucl. Acids Res. 41 (D1): D590-D596.
- Yilmaz P, Parfrey LW, Yarza P, Gerken J, Pruesse E, Quast C, Schweer T, Peplies J, Ludwig W, Glöckner FO (2014) The SILVA and "All-species Living Tree Project (LTP)" taxonomic frameworks. Nucl. Acids Res. 42:D643-D648
- Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. 2016. DADA2: High-resolution sample inference from Illumina amplicon data. Nat Methods 13:581–583. doi:10.1038/nmeth.3869

##### UNITE
- Nilsson, R. H. et al. (2019) The UNITE database for molecular identification of fungi: handling dark taxa and parallel taxonomic classifications. Nucleic Acids Res 47, D259–D264.
