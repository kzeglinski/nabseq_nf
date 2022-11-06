#!/usr/bin/env nextflow

/*
 * NAb-seq nextflow pipeline for the analysis of Oxford Nanopore sequencing data from
 * hybridomas and single B cells 
 *
 * If you use NAb-seq, please cite https://www.tandfonline.com/doi/full/10.1080/19420862.2022.2106621
 *
 */

version = "v0.2"

// TODO: help message
log.info """
##    ##    ###    ########           ######  ########  #######  
###   ##   ## ##   ##     ##         ##    ## ##       ##     ## 
####  ##  ##   ##  ##     ##         ##       ##       ##     ## 
## ## ## ##     ## ########  #######  ######  ######   ##     ## 
##  #### ######### ##     ##               ## ##       ##  ## ## 
##   ### ##     ## ##     ##         ##    ## ##       ##    ##  
##    ## ##     ## ########           ######  ########  ##### ## 

================================================================
Usage: nextflow run ./nabseq_nf/main.nf --fastq_dir [input path] --organism [organism name] --sample_sheet [sample sheet]
--help                : prints this help message

Required arguments:
--out_dir             : where the output files will be written to (default: "$projectDir/results)
--fastq_dir           : where the input fastq files are located
--organism            : name of the organism whose reference sequences to use (default: "rat")
--sample_sheet        : location of the .csv sample sheet (format: barcode01,sample_x)

Optional (only needed for advanced users)
--num_consensus       : maximum number of consensus sequences to generate for each sample (default: 999)
--igblast_databases   : location of the igblast databases (default: "$projectDir/references/igblast/")
--reference_sequences : location of the reference sequences (default: "$projectDir/references/reference_sequences/")
--trim_3p             : pattern for cutadapt to trim off the 3' end (default: "A{20}N{90}")
--trim_5p             : pattern for cutadapt to trim off the 5' end (default: "N{90}T{20}")
--medaka_model        : model to use for medaka (depends on your basecalling model, default:"r941_min_sup_g507")
--report_title        : title to use for naming the report (default: "NAb-seq report")

"""

// TODO: parameter validation
// validate that these files/directories exist
fastq_dir = params.fastq_dir
sample_sheet = params.sample_sheet
reference_dir = params.reference_sequences
igblast_databases = params.igblast_databases
igdata_dir="${igblast_databases}/igdata/"
igblastdb_dir="${igblast_databases}/databases/"

// validate that these are strings 
organism = params.organism // check that there are references for this organism
trim_3p = params.trim_3p
trim_5p = params.trim_5p
num_consensus = params.num_consensus
medaka_model = params.medaka_model
report_title = params.report_title

// TODO: introductory message that prints out parameters
log.info """
##    ##    ###    ########           ######  ########  #######  
###   ##   ## ##   ##     ##         ##    ## ##       ##     ## 
####  ##  ##   ##  ##     ##         ##       ##       ##     ## 
## ## ## ##     ## ########  #######  ######  ######   ##     ## 
##  #### ######### ##     ##               ## ##       ##  ## ## 
##   ### ##     ## ##     ##         ##    ## ##       ##    ##  
##    ## ##     ## ########           ######  ########  ##### ## 

================================================================
read directory              : ${params.fastq_dir}
number of consensus seqs    : ${params.num_consensus}
organism                    : ${params.organism}
sample_sheet                : ${params.sample_sheet}
medaka_model                : ${params.medaka_model}
output directory            : ${params.out_dir}
report_title                : ${params.report_title}

"""

/*
 * Bring in modules
 */
include { parse_sample_sheet } from './subworkflows/file_import'
include { nanocomp } from './modules/local/nanocomp' 
include { select_ab_reads } from './subworkflows/select_ab_reads'
include { read_trimming } from './subworkflows/read_trimming'
include { fastq_to_fasta } from './modules/local/fastq_to_fasta'
include { annotation_grouping_pre_consensus } from './subworkflows/annotation_grouping_pre_consensus'
include { take_consensus } from './subworkflows/consensus'
include { annotation_grouping_post_consensus } from './subworkflows/annotation_grouping_post_consensus'
include { report } from './subworkflows/report'


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
    
    // trim the ends (by default polyA tail, user can specify if they have primers)
    trimmed_reads = read_trimming(ab_reads, trim_3p, trim_5p).trimmed_reads_renamed    

    // annotation and grouping
    trimmed_read_fasta = fastq_to_fasta(trimmed_reads) //IgBLAST needs fasta not fastq
    annotation_grouping_pre_consensus(
        trimmed_read_fasta, 
        organism, 
        igblast_databases, 
        igdata_dir, 
        igblastdb_dir,
        num_consensus)
    pre_consensus_table = annotation_grouping_pre_consensus.out.pre_consensus_table
    
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
    .map{ it -> [it[0].replaceAll('_[^_]*_(.*)', "")] + it} // add sample ID to list
    .combine(trimmed_reads, by: 0) // add trimmed reads to list
    .set{input_for_consensus}

     // actually make the consensus 
    consensus_sequences = take_consensus(input_for_consensus, medaka_model)
    
    // annotation & grouping/analysis of consensus sequences
    // as input for this, we need the pre-consensus grouped table, as well as the consensus sequences
    consensus_sequences
    .map{ it -> [it[0].replaceAll('_[^_]*_(.*)', "")] + it }
    .groupTuple() 
    .set{final_annotation_files}

    annotation_grouping_post_consensus(
        final_annotation_files,
        organism,
        igblast_databases,
        igdata_dir,
        igblastdb_dir)
    
    productive_only_annotation = annotation_grouping_post_consensus.out.only_productive_sequences.collect()
    full_annotation = annotation_grouping_post_consensus.out.all_annotated_sequences.collect()

    // make report
    report_template = Channel.fromPath(params.report_template)
    css_file = Channel.fromPath(params.css_file)
    logo = Channel.fromPath(params.logo_file)
    
    nanocomp.out.png
    .flatten()
    .concat(nanocomp.out.txt)
    .collect()
    .set{nanocomp_output}
    
    report(
        version,
        organism,
        sample_sheet,
        report_title,
        report_template,
        productive_only_annotation,
        full_annotation,
        nanocomp_output,
        css_file,
        logo
    )

}