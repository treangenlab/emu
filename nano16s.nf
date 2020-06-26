#!/usr/bin/env nextflow

params.in = "/scratch/Maria/16s_full_length_mock/nano16s_nextflow/SRR8029999.fastq"
// params.ref = "/scratch/Maria/16s_full_length_mock/rrnDB-5.6_16S_rRNA_modification.fasta"
params.ref = "/scratch/Maria/silva_db.fasta"
params.dict = "/scratch/Maria/16s_full_length_mock/nano16s_nextflow/templates/silva_tax_dict"
params.outputDir = "/scratch/Maria/16s_full_length_mock/nano16s_nextflow/results"
query_ch = Channel.fromPath(params.in)

process porechop{
    
    publishDir './results'
    
    input:
    file query from query_ch

    output:
    file 'porechop.fastq' into output
    

    script:
        """
	porechop -t 40 -i $query -o porechop.fastq """
} 

process length_selection{

    publishDir './results'

    input:
    file x from output

    output:
    file('demux_len_1200.fasta') into selection

    script:
    """
    awk 'BEGIN {FS="\t"; OFS="\\n"} {header = \$0 ; getline seq ; getline qheader ; getline qseq ; if (length(seq) >= 1200) {print header, seq, qheard, qseq}}' < $x > demux_len_1200.fastq
    sed '/^@/!d;s//>/;N' demux_len_1200.fastq > demux_len_1200.fasta 
    rm demux_len_1200.fastq
    """
}

process minimap{
    publishDir './results'
    
    input:
    file y from selection

    output:
    file("demux_len_1200_silva_map-ont") into mapping
    
    script:
    """
    minimap2 -x map-ont -z 70 -t 40 $params.ref $y -o demux_len_1200_silva_map-ont
    """
}

process parse_minimap{
    publishDir './results'

    input:
    file(alignment) from mapping

    output:
    file("tax_sep_silva_sample.tsv") into parsing

    script:
    template "parse_minimap2_silva.py"
}


process get_abundance{
    publishDir './results'

    input:
    file(reads) from parsing
    
    output:
    tuple file("abundance_silva_genus.tsv"), file("abundance_silva_sp.tsv")

    script:
    template "get_abundance.py"
}
