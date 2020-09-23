# Nano16s

full-length 16s classification software


## Expectation Maximization Classification

EM Classification algorithm for long-read 16s sequences.

output: .tsv file containing lineages and relative abundances for input file in `results/` dir

```bash
python EM_classification.py <input_file [fasta, fastq, sam]>
```

```bash
python EM_classification.py example/basic4.fa
```
