// this is the pre consensus annotation and grouping
// convert fastq to fasta, then run IgBLAST
// process this in R

include { igblast } from '../modules/local/igblast'

process cat_all_abs{
    tag "$prefix"
    label 'process_low'
    publishDir "${params.out_dir}/consensus_annotation", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(meta), path(consensus_sequences)

    output:
    tuple val(meta), path("*_all_consensus_sequences.fasta")

    script:
    // allow for a bunch of metadata (although the first element should be sample name)
    if(meta instanceof Collection) {
        prefix = meta[0]
    } else {
        prefix = meta
    }

    """
    cat $consensus_sequences > "${prefix}_all_consensus_sequences.fasta"
    """
}

process post_consensus_annotation {
    tag "$prefix"
    label 'process_low'
    publishDir "${params.out_dir}/consensus_annotation", mode: 'copy', pattern: "*.tsv"
    container "library://kzeglinski/nabseq/nabseq-report:v0.0.3"

    input:
    tuple val(meta), path(igblast_output)
    path(sample_sheet)

    output:
    //tuple val(meta), path('*_full_consensus_annotation.tsv'), emit: full_annotation
    tuple val(meta), path('*_productive_only_consensus_annotation.tsv'), emit: prod_annotation
    tuple val(meta), path('*_flag_data.tsv'), emit: flag_data

    script:
    // allow for a bunch of metadata (although the first element should be sample name)
    if(meta instanceof Collection) {
        prefix = meta[0]
    } else {
        prefix = meta
    }

    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(tidyr)
    library(stringr)
    library(vroom)
    library(Biostrings)

    # read in the igblast output
    igblast_results <- vroom("${igblast_output}",
        col_select = c(
            sequence_id, sequence, locus, productive, complete_vdj, v_call,
            d_call, j_call, c_call, cdr3_aa, v_sequence_start, c_sequence_end))

    if(nrow(igblast_results) == 0){
        stop("No IgBLAST results found")
    }

    # separate out the counts from the name
    igblast_results %>%
        extract(sequence_id, into = c("sequence_id", "count"), "(.*)_([^_]+)\$") %>%
        mutate(group_id = gsub('_count', '', sequence_id), .before = count) %>%
        arrange(desc(count)) %>%
        select(-sequence_id) -> igblast_results_with_counts

    # trim to just the V(D)JC region unless otherwise specified
    if(as.logical("${params.trim_to_antibody}")){
        igblast_results_with_counts %>%
            mutate(sequence_untrimmed = sequence) %>%
            mutate(sequence = str_sub(
                sequence,
                start = v_sequence_start,
                end = ifelse(!is.na(c_sequence_end), c_sequence_end, -1))) -> igblast_results_with_counts
    }

    # write out a copy of this table
    vroom_write(igblast_results_with_counts, "${prefix}_full_consensus_annotation.tsv")

    # just productive ones
    igblast_results_with_counts %>%
        filter(productive == TRUE) %>%
        select(-sequence_untrimmed) %>%
        group_by(v_call, j_call, cdr3_aa, c_call) %>%
        summarise(
            group_id = paste0(group_id, collapse = ", "),
            count = sum(as.numeric(count)),
            locus = locus[1],
            d_call = d_call[1],
            sequence = sequence[which.max(nchar(sequence))]) -> igblast_results_prod_only

    vroom_write(igblast_results_prod_only, "${prefix}_productive_only_consensus_annotation.tsv")

    # prepare flag data for report
    flag_data <- data.frame(
        sample_id = str_remove(igblast_results_prod_only[["group_id"]][1], "_.*"),
        chain_status = "",
        sp2_0_status = "normal"
    )

    # cover the case where we have no productive sequences
    if(is.na(flag_data[["sample_id"]][1])){
        flag_data[["sample_id"]] <- str_remove(igblast_results_with_counts[["group_id"]][1], "_.*")
    }

    # check for multiple H/L chains
    igblast_results_prod_only %>%
        filter(!is.na(cdr3_aa)) %>% # dont flag duplicates with no CDR3
        ungroup() %>%
        group_by(cdr3_aa) %>%
        dplyr::slice_max(order_by = count, n = 1) -> dup_cdr3_removed # don't flag duplicates with the same CDR3

    if(length(unique(dup_cdr3_removed[["locus"]])) < 2){
        flag_data[["chain_status"]] <- paste(flag_data[["chain_status"]], "Missing H/L chain")
    }

    if(any(duplicated(dup_cdr3_removed[["locus"]]))){
        flag_data[["chain_status"]] <- paste(flag_data[["chain_status"]], " Multiple H/L chains")
    }

    if(flag_data[["chain_status"]] == ""){
        flag_data[["chain_status"]] <- "normal"
    }

    # remove whitespace
    flag_data[["chain_status"]] <- str_trim(flag_data[["chain_status"]])

    # align to the aberrant sp2/0 light chain
    sp2_0 <- DNAString("ATGGAGACAGACACACTCCTGTTATGGGTACTGCTGCTCTGGGTTCCAGGTTCCACTGGTGACATTGTGCTGACACAGTCTCCTGCTTCCTTAGCTGTATCTCTGGGGCAGAGGGCCACCATCTCATACAGGGCCAGCAAAAGTGTCAGTACATCTGGCTATAGTTATATGCACTGGAACCAACAGAAACCAGGACAGCCACCCAGACTCCTCATCTATCTTGTATCCAACCTAGAATCTGGGGTCCCTGCCAGGTTCAGTGGCAGTGGGTCTGGGACAGACTTCACCCTCAACATCCATCCTGTGGAGGAGGAGGATGCTGCAACCTATTACTGTCAGCACATTAGGGAGCTTACACGTTCGGAGGGGGGACCAAGCTGGAAATAAAACGGGCTGATGCTGCACCAACTGTATCCATCTTCCCACCATCCAGTGAGCAGTTAACATCTGGAGGTGCCTCAGTCGTGTGCTTCTTGAACAACTTCTACCCCAAAGACATCAATGTCAAGTGGAAGATTGATGGCAGTGAACGACAAAATGGCGTCCTGAACAGTTGGACTGATCAGGACAGCAAAGACAGCACCTACAGCATGAGCAGCACCCTCACGTTGACCAAGGACGAGTATGAACGACATAACAGCTATACCTGTGAGGCCACTCACAAGACATCAACTTCACCCATTGTCAAGAGCTTCAACAGGAATGAGTGTTAGAGACAAAGGTCCTGAGACGCCACCACCAGCTCCCCAGCTCCATCCTATCTTCCCTTCTAAGGTCTTGGAGGCTTCCCCACAAGCGACCTACCACTGTTGCGGTGCTCCAAACCTCCTCCCCACCTCCTTCTCCTCCTCCTCCCTTTCCTTGGCTTTTATCATGCTAATATTTGCAGAAAATATTCAATAAAGTGAGTCTTTGCACTTGAAAAAAAAAAAAA")

    align_these <- DNAStringSet(igblast_results_prod_only[["sequence"]])
    sp2_0_identities <- pid(pairwiseAlignment(align_these, sp2_0))

    if(any(sp2_0_identities > 90)){
        flag_data[["sp2_0_status"]] <- "Aberrant sp2/0 light chain detected"
    }

    # prepare the final flag
    if(flag_data[["sp2_0_status"]] == "normal" & flag_data[["chain_status"]] == "normal"){
        flag_data[["flag"]] <- "normal"
    } else if(flag_data[["sp2_0_status"]] == "normal"){
        flag_data[["flag"]] <- flag_data[["chain_status"]]
    } else if(flag_data[["chain_status"]] == "normal") {
        flag_data[["flag"]] <- flag_data[["sp2_0_status"]]
    } else {
        flag_data[["flag"]] <- paste(
            flag_data[["chain_status"]], flag_data[["sp2_0_status"]], sep = ",")
    }

    # write out the flag data
    vroom_write(flag_data, "${prefix}_flag_data.tsv")
    """
}

workflow annotation_grouping_post_consensus {
    take:
        final_annotation_files
        igblast_databases
        sample_sheet

    main:
        // cat all consensus sequences from the same sample
        concatenated_sequences = cat_all_abs(final_annotation_files)

        // set up environment variables
        igdata_dir="${igblast_databases}/igdata/"
        igblastdb_dir="${igblast_databases}/databases/"

        // annotate reads using igblast
        igblast_tsv = igblast(
            concatenated_sequences,
            igblast_databases,
            igdata_dir,
            igblastdb_dir,
            "post").airr_table

        // process this in R
        post_consensus_annotation(igblast_tsv, sample_sheet)

        // prepare output files
        annotated_sequences = post_consensus_annotation.out.prod_annotation
        flag_data = post_consensus_annotation.out.flag_data


    emit:
        annotated_sequences
        flag_data

}