// this is the pre consensus annotation and grouping
// convert fastq to fasta, then run IgBLAST
// process this in R 

include { igblast } from '../modules/local/igblast' 

workflow annotation_grouping_pre_consensus {
    take:
        trimmed_fasta
        organism
        igblast_databases
        igdata_dir
        igblastdb_dir
        
    
    main:
        // annotate reads using igblast
        igblast_tsv = igblast(trimmed_fasta, organism, igblast_databases, igdata_dir, igblastdb_dir).airr_table

        // process this in R

        

    //emit:
    //    ab_reads

}