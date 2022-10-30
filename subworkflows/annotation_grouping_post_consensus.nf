// this is the pre consensus annotation and grouping
// convert fastq to fasta, then run IgBLAST
// process this in R 

include { igblast } from '../modules/local/igblast' 

process cat_all_abs{
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.out_dir}/consensus_annotation", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(sample_id), val(sequence_ids), path(consensus_sequences), path(pre_consensus_table)

    output:
    tuple val(sample_id), path("*_all_consensus_sequences.fasta"), path(pre_consensus_table)

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
    tuple val(meta), path(igblast_output), path(consensus_sequence), path(pre_consensus_table)

    output:
    tuple val(meta), path('*_full_consensus_annotation.tsv'), emit: full_annotation
    tuple val(meta), path('*_productive_only_consensus_annotation.tsv'), emit: prod_annotation

    script:
    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    
    # read in the igblast output
    igblast_results <- read.delim(file = "${igblast_output}", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

    # also read in the pre-consensus table (since this has the counts)
    pre_consensus_igblast <- read.delim(file = "${pre_consensus_table}", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

    pre_consensus_igblast %>%
        select(c(count, group_id)) %>%
        right_join(igblast_results, by = c("group_id" = "sequence_id")) -> igblast_results_with_counts

    # write out a copy of this table 
    write_tsv(igblast_results_with_counts, "${meta}_full_consensus_annotation.tsv")

    # just productive ones
    # remove those that have less than 3 counts (can't make a consensus with 1 or 2 reads)
    igblast_results_with_counts %>%
        filter(productive == TRUE) %>%
        select(c(group_id, productive, v_call, d_call, j_call, cdr3_aa, v_identity, c_call, count))-> igblast_results_prod_only
    
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
        concatenated_sequences.map{it -> it.take(2)}.set{igblast_input} // don't want to give the .tsv to igblast

        // annotate reads using igblast
        igblast_tsv = igblast(
            igblast_input, 
            organism, 
            igblast_databases, 
            igdata_dir, 
            igblastdb_dir,
            "post").airr_table
        
        // add on the pre-consensus table so we can do the final analysis in R
        igblast_tsv
        .combine(concatenated_sequences, by: 0)
        .set{for_final_processing}

        // process this in R
        post_consensus_annotation(for_final_processing)
        all_annotated_sequences = post_consensus_annotation.out.full_annotation
        only_productive_sequences = post_consensus_annotation.out.prod_annotation

    emit:
        all_annotated_sequences
        only_productive_sequences

}