# Nano16s

full-length 16s classification software


## Expectation Maximization Classification

**dependencies:**




EM Classification algorithm for long-read 16s sequences. Uses NCBI 16s refseq database.

Input file must be in one of the following formats: fasta, fastq, sam.

* fasta/q: performs BWA parameterized for Oxford Nanopore Technologies long read data, then performs classification algorithm on output sam file.
* sam: file contains all BWA alignments for input reads compared to NCBI 16s refseq database. Performs classification algorithm on input sam file.

output: .tsv file containing lineages and relative abundances for input file in `results/` dir

```bash
python EM_classification.py <input_file [fasta, fastq, sam]>
```

```bash
python EM_classification.py example/basic4.fa
```
