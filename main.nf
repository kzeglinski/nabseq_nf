#!/usr/bin/env nextflow

// TODO: intro/copyright notice

// TODO: help message

// TODO: parameter validation
// TODO: check these files/folders exist
fastq_dir = params.fastq_dir
sample_sheet = params.sample_sheet
reference_dir = params.reference_sequences
igblast_databases = params.igblast_databases
igdata_dir="${igblast_databases}/igdata/"
igblastdb_dir="${igblast_databases}/databases/"

// TODO: check these are strings
organism = params.organism
trim_3p = params.trim_3p
trim_5p = params.trim_5p
num_consensus = params.num_consensus

// TODO: introductory message that prints out parameters



/*
 * Bring in modules
 */
include { parse_sample_sheet } from './subworkflows/file_import'
include { nanocomp } from './modules/local/nanocomp' 
include { select_ab_reads } from './subworkflows/select_ab_reads'
include { cutadapt } from './modules/local/cutadapt'
include { fastq_to_fasta } from './modules/local/fastq_to_fasta'
include { annotation_grouping_pre_consensus } from './subworkflows/annotation_grouping_pre_consensus'
include { take_consensus } from './subworkflows/consensus'

/*
 * Define helper functions
 */

// the last part of this is inspired by https://github.com/stevekm/nextflow-demos/blob/master/join-pairs/main.nf
def reshape_channel_for_consensus (input_channel, file_ending) {
    input_channel
    .map { it -> [it.name -  ~/\.txt/, it] } // get the file names to get the ID from
    .collect(flat: false){ it[0] } // subset out those filenames
    .collectNested { it.minus(file_ending) } // remove the file name so it's just the ID
    .combine(input_channel) // combine the lists so we can pair 
    .flatten().toList() // honestly idk why but this combination gets it into a form the pairs() function works on
    .flatMap { pairs(it) } // apply the pairs function
    .filter { items -> // only keep combinations where the ID matches the filename
            def id_1 = items[0].toString().replaceAll('.*/{1}', '').replaceAll(file_ending, "") //strips path
            def id_2 = items[1].toString().replaceAll('.*/{1}', '').replaceAll(file_ending, "")
            id_1 == id_2
            }
    .filter { items -> // remove the pairs that are like [ID, ID] instead of [ID, file] or [file, ID]
            def id_1 = items[0]
            def id_2 = items[1]
            id_1 != id_2
            }
    .filter { items -> // currently we have [ID, file] and [file, ID] but we only need one of them
            def id_1 = items[0].toString()
            !id_1.endsWith(".txt") // remove [file, ID] so we're just left with [ID, file]
            }
}

// from https://stackoverflow.com/questions/27868047/groovy-way-to-pair-each-element-in-a-list
// this is needed to reshape the input channel for consensus calling
def pairs(def elements) {
    return elements.tail().collect { [elements.head(), it] } + (elements.size() > 1 ? pairs(elements.tail()) : [])
}

/*
 * Run the workflow
 */

    
workflow{
    // process the sample sheet, then concat all fastq files in
    // the barcode directories & name based on the sample (not barcode)
    concatenated_files = parse_sample_sheet(fastq_dir, sample_sheet)

    // since the channel from the previous step includes sample names, just get fastq
    // then collect and send to nanocomp for QC
    nanocomp_input = concatenated_files.flatten().filter(~/.*fastq.gz$/).collect() 
    nanocomp(nanocomp_input)

    // select antibody reads (w/ minimap2 alignment)
    all_ab_reference_file = file("${reference_dir}/${organism}_all_imgt_refs.fasta")
    if( !all_ab_reference_file.exists() ) exit 1, "The reference file ${all_ab_reference_file} can't be found."
    ab_reads = select_ab_reads(concatenated_files, all_ab_reference_file)
    
    // trim the 3' and 5' ends (by default polyA tail, user can specify if they have primers)
    trimmed_reads = cutadapt(ab_reads, trim_3p, trim_5p).reads    

    // annotation and grouping
    trimmed_read_fasta = fastq_to_fasta(trimmed_reads) //IgBLAST needs fasta not fastq
    annotation_grouping_pre_consensus(
        trimmed_read_fasta, 
        organism, 
        igblast_databases, 
        igdata_dir, 
        igblastdb_dir,
        num_consensus)
    
    // the .collect(flat: false) doesn't do anything to the output. It just means that the next step (reshaping
    // the channel) won't start until annotation_grouping_pre_consensus is finished running on all samples
    // otherwise you end up with a shitload of duplicates in the next step
    read_names_for_consensus = annotation_grouping_pre_consensus.out.read_names_for_consensus.collect(flat: false)
    starting_points_for_consensus = annotation_grouping_pre_consensus.out.starting_points_for_consensus.collect(flat: false)

    // consensus calling 
    // getting the channels into the right shape
    // this section is so clunky and messy and awful but i don't know enough to make it better
    reshaped_read_names = reshape_channel_for_consensus(read_names_for_consensus, "_read_names.txt")
    reshaped_starting_points = reshape_channel_for_consensus(starting_points_for_consensus, "_starting_point_name.txt")

    reshaped_read_names
    .join(reshaped_starting_points)
    .map{ it -> [it[0].replaceAll('_[^_]+$', "")] + it} // add sample ID to list
    .combine(trimmed_reads, by: 0) // add trimmed reads to list
    .set{input_for_consensus}

    // actually make the consensus 
    consensus_sequences = take_consensus(input_for_consensus)
    
    // annotation & grouping/analysis of consensus sequences

    // make report

    // finished! print some kind of message to user
    // also email if they want?
}