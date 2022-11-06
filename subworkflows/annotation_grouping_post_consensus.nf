// this is the pre consensus annotation and grouping
// convert fastq to fasta, then run IgBLAST
// process this in R 

include { igblast } from '../modules/local/igblast' 

process cat_all_abs{
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.out_dir}/consensus_annotation", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(sample_id), val(sequence_ids), path(consensus_sequences)

    output:
    tuple val(sample_id), path("*_all_consensus_sequences.fasta")

    script:

    """
    cat $consensus_sequences > "${sample_id}_all_consensus_sequences.fasta"
    """
}

process post_consensus_annotation {
    tag "$meta"
    label 'process_low'
    publishDir "${params.out_dir}/consensus_annotation", mode: 'copy', pattern: "*.tsv"

    conda (params.enable_conda ? 'r::r-tidyverse=1.2.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse%3A1.2.1' :
        'quay.io/biocontainers/r-tidyverse:1.2.1' }"

    input:
    tuple val(meta), path(igblast_output)

    output:
    path('*_full_consensus_annotation.tsv'), emit: full_annotation
    path('*_productive_only_consensus_annotation.tsv'), emit: prod_annotation

    script:
    
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    
    # read in the igblast output
    igblast_results <- read.delim(file = "${igblast_output}", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

    # separate out the counts from the name
    igblast_results %>% 
        extract(sequence_id, into = c("sequence_id", "count"), "(.*)_([^_]+)\$") %>%
        mutate(group_id = gsub('_count', '', sequence_id), .before = count) %>%
        arrange(desc(count)) %>%
        select(-sequence_id) -> igblast_results_with_counts

    # write out a copy of this table 
    write_tsv(igblast_results_with_counts, "${meta}_full_consensus_annotation.tsv")

    # just productive ones
    # remove those that have less than 3 counts (can't make a consensus with 1 or 2 reads)
    igblast_results_with_counts %>%
        filter(productive == TRUE) %>%
        select(c(group_id, count, productive, v_call, d_call, j_call, cdr3_aa, v_identity, c_call)) %>%
        ungroup() %>%
        group_by(v_call, d_call, j_call, cdr3_aa, v_identity, c_call) %>%
        summarise(group_id = paste0(group_id, collapse = ", "), count = sum(as.numeric(count))) -> igblast_results_prod_only
    
    write_tsv(igblast_results_prod_only, "${meta}_productive_only_consensus_annotation.tsv")
        
    """
}

workflow annotation_grouping_post_consensus {
    take:
        final_annotation_files
        organism
        igblast_databases
        igdata_dir
        igblastdb_dir
        
    main:
        // cat all consensus sequences from the same sample
        concatenated_sequences = cat_all_abs(final_annotation_files)

        // annotate reads using igblast
        igblast_tsv = igblast(
            concatenated_sequences, 
            organism, 
            igblast_databases, 
            igdata_dir, 
            igblastdb_dir,
            "post").airr_table
        
        // process this in R
        post_consensus_annotation(igblast_tsv)
        all_annotated_sequences = post_consensus_annotation.out.full_annotation
        only_productive_sequences = post_consensus_annotation.out.prod_annotation

    emit:
        all_annotated_sequences
        only_productive_sequences

}