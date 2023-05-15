// this is the pre consensus annotation and grouping
// convert fastq to fasta, then run IgBLAST
// process this in R 

include { igblast } from '../modules/local/igblast' 

process pre_consensus_groupings {
    tag "$prefix"
    label 'process_low'
    publishDir "${params.out_dir}/pre_consensus_grouping", mode: 'copy', pattern: "*.tsv"

    conda (params.enable_conda ? 'r::r-tidyverse=1.2.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse%3A1.2.1' :
        'quay.io/biocontainers/r-tidyverse:1.2.1' }"

    input:
    tuple val(meta), path(igblast_output)
    val num_consensus

    output:
    tuple val(meta), path('*_pre_consensus_grouped_table.tsv'), emit: grouped_table
    path("*_starting_point_name.txt"), emit: starting_points
    path("*_read_names.txt"), emit: read_names

    script:
    // allow for a bunch of metadata (although the first element should be sample name)
    if(meta instanceof Collection) {
        prefix = meta[0]
    } else {
        prefix = meta
    }

    """
    #!/usr/bin/env Rscript

    library(tidyverse)
    # there should only be one .tsv file per folder, so we can
    # read it in like so:
    igblast_results <- read.delim(file = list.files(pattern='.tsv')[1], sep = "\t", header = TRUE, stringsAsFactors = FALSE)

    # add n column to count up how many reads fall into each
    # then group the reads by V(D)JC
    igblast_results %>%
    select(c(v_call, d_call, j_call, c_call, sequence_id)) %>%
    mutate(n = rep(1, nrow(igblast_results))) %>%
    group_by(v_call, d_call, j_call, c_call) %>%
        summarise(count = sum(n), .groups = "keep", 
            reads = paste(sequence_id, collapse = "_")) -> igblast_results_grouped

    # separate out the heavy and light chains
    igblast_results_grouped_H_full <- filter(igblast_results_grouped, str_detect(v_call, "H"))
    igblast_results_grouped_L_full <- filter(igblast_results_grouped,
                                             str_detect(v_call, "H", negate = TRUE))

    # sort them in descending order by read count
    igblast_results_grouped_H_full <- arrange(igblast_results_grouped_H_full, desc(count))
    igblast_results_grouped_L_full <- arrange(igblast_results_grouped_L_full, desc(count))

    # give each clones H and L chain a unique name in the form of H1, H2, L1, L2 etc where 1 is the most abundant, 2 is the second most abundant etc
    igblast_results_grouped_H_full[, "group_id"] <- paste0("$prefix", "_H", seq_len(nrow(igblast_results_grouped_H_full)), "_count_", as.character(unlist(as.vector(igblast_results_grouped_H_full[, "count"]))))
    igblast_results_grouped_L_full[, "group_id"] <- paste0("$prefix", "_L", seq_len(nrow(igblast_results_grouped_L_full)), "_count_", as.character(unlist(as.vector(igblast_results_grouped_L_full[, "count"]))))

    # write out a copy of this table 
    igblast_results_grouped_full <- bind_rows(igblast_results_grouped_H_full, igblast_results_grouped_L_full)
    write_tsv(igblast_results_grouped_full, "${prefix}_pre_consensus_grouped_table.tsv")

    # prepare for consensus calling 
    # remove those that have less than 3 counts (can't make a consensus with 1 or 2 reads)
    igblast_results_grouped_H_full %>%
        filter(count >= 3) -> igblast_results_grouped_H

    igblast_results_grouped_L_full %>%
        filter(count >= 3) -> igblast_results_grouped_L

    # select the top n groups, based on the parameter set
    # default value is 999 which should keep all 
    if ($num_consensus < nrow(igblast_results_grouped_H)) {
        igblast_results_grouped_H <- igblast_results_grouped_H[1:$num_consensus,]
    }

    if ($num_consensus < nrow(igblast_results_grouped_L)) {
        igblast_results_grouped_L <- igblast_results_grouped_L[1:$num_consensus,]
    }

    # write out the read names for each clone into a .txt file
    for (i in seq_along(unlist(as.vector(igblast_results_grouped_H[, "group_id"])))) {
        this_group_id <- unlist(as.vector(igblast_results_grouped_H[, "group_id"]))[i]
        this_reads <- unlist(str_split(unlist(as.vector(igblast_results_grouped_H[, "reads"]))[i], "_"))
        this_file <- paste0(this_group_id, "_read_names.txt")
        writeLines(this_reads, this_file)
    }

    for (i in seq_along(unlist(as.vector(igblast_results_grouped_L[, "group_id"])))) {
        this_group_id <- unlist(as.vector(igblast_results_grouped_L[, "group_id"]))[i]
        this_reads <- unlist(str_split(unlist(as.vector(igblast_results_grouped_L[, "reads"]))[i], "_"))
        this_file <- paste0(this_group_id, "_read_names.txt")
        writeLines(this_reads, this_file)
    }

    # choose the starting copy 
    # make a long format data where each read is a row (but still keep track of which group each read belongs to) 
    igblast_results_grouped_long <- bind_rows(
        separate_rows(igblast_results_grouped_L, reads, sep = "_"),
        separate_rows(igblast_results_grouped_H, reads, sep = "_"))

    # get the longest read in each group (since it's more likely to have a full-length C region)
    # find read length from the nucleotide sequence
    # first just keep only essential columns
    igblast_results_grouped_long %>%
        select(c(reads, group_id)) -> igblast_results_grouped_long

    # get nucleotide sequence from original table
    # join them, remove those rows that don't belong to a group and determine read length
    # finally, just choose the longest reads for each group
    # EDIT: [21/12/22] we actually need to choose the longest read that also has a complete VDJ
    # otherwise we risk outputting a truncated consensus sequence
    igblast_results %>%
        select(c(sequence_id, sequence, complete_vdj)) %>%
        left_join(igblast_results_grouped_long, by = c("sequence_id" = "reads")) %>%
        filter(!is.na(group_id)) %>%
        filter(complete_vdj == TRUE) %>%
        mutate(read_length = nchar(sequence)) %>%
        select(-c(sequence)) %>%
        group_by(group_id) %>% 
        slice(which.max(read_length)) -> igblast_results_grouped_longest_complete_reads
    
    # but sometimes we might not have any complete VDJ for a group
    # in that case, just choose the longest read
    igblast_results %>%
        select(c(sequence_id, sequence)) %>%
        left_join(igblast_results_grouped_long, by = c("sequence_id" = "reads")) %>%
        filter(!is.na(group_id)) %>%
        mutate(read_length = nchar(sequence)) %>%
        select(-c(sequence)) %>%
        group_by(group_id) %>% 
        slice(which.max(read_length)) -> igblast_results_grouped_longest_reads

    # remove the rows of igblast_results_grouped_longest_reads that appear in igblast_results_grouped_longest_complete_reads
    igblast_results_grouped_longest_reads %>%
        anti_join(igblast_results_grouped_longest_complete_reads, by = c("sequence_id" = "sequence_id")) -> igblast_results_grouped_longest_reads_no_dupes

    # then can combine the two tables to come up with our starting point master list
    igblast_results_grouped_longest_complete_reads %>%
        bind_rows(igblast_results_grouped_longest_reads_no_dupes) -> starting_point_reads
    
    # write out these longest reads as the starting copies
    for (i in seq_along(unlist(as.vector(starting_point_reads[, "group_id"])))) {
        this_group_id <- unname(unlist(as.vector(starting_point_reads[, "group_id"]))[i])
        this_read <- unname(unlist(as.vector(starting_point_reads[, "sequence_id"]))[i])
        this_file <- paste0(this_group_id, "_starting_point_name.txt")
        writeLines(this_read, this_file)
    }
        
    """
}

workflow annotation_grouping_pre_consensus {
    take:
        trimmed_fasta
        igblast_databases
        num_consensus
        
    main:
        // set up environment variables
        igdata_dir="${igblast_databases}/igdata/"
        igblastdb_dir="${igblast_databases}/databases/"

        // annotate reads using igblast
        igblast_tsv = igblast(
            trimmed_fasta, 
            igblast_databases, 
            igdata_dir, 
            igblastdb_dir,
            "pre").airr_table

        // process this in R
        pre_consensus_groupings(igblast_tsv, num_consensus)
        read_names_for_consensus = pre_consensus_groupings.out.read_names
        starting_points_for_consensus = pre_consensus_groupings.out.starting_points
        pre_consensus_table = pre_consensus_groupings.out.grouped_table

    emit:
        read_names_for_consensus
        starting_points_for_consensus
        pre_consensus_table

}
