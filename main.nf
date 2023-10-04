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
if(params.help == true){
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
--sample_sheet        : location of the .csv sample sheet (format: barcode01,sample_x,rat,1)

Optional (only needed for advanced users)
--barcode_dirs        : whether the input fastq files are located within folders named barcode01 etc (default: false)
--num_consensus       : maximum number of consensus sequences to generate for each sample (default: 999)
--igblast_databases   : location of the igblast databases (default: "$projectDir/references/igblast/")
--reference_sequences : location of the reference sequences (default: "$projectDir/references/reference_sequences/")
--trim_3p             : pattern for cutadapt to trim off the 3' end (default: "A{20}N{90}")
--trim_5p             : pattern for cutadapt to trim off the 5' end (default: "N{90}T{20}")
--medaka_model        : model to use for medaka (depends on your basecalling model, default:"r941_min_sup_g507")
--report_title        : title to use for naming the report (default: "NAb-seq report")

"""
System.exit(0)
}
// TODO: parameter validation
// validate that these files/directories exist
fastq_dir = params.fastq_dir
sample_sheet = params.sample_sheet
reference_dir = params.reference_sequences
igblast_databases = params.igblast_databases

// validate that these are strings 
trim_3p = params.trim_3p
trim_5p = params.trim_5p
num_consensus = params.num_consensus
medaka_model = params.medaka_model
report_title = params.report_title
report_template = params.report_template
css_file = params.css_file
logo = params.logo_file

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
include { versions } from './subworkflows/versions'
include { validate_params } from './subworkflows/validate_params'
/*
 * Run the workflow
 */
workflow{

    pathsToValidate = [fastq_dir, sample_sheet, reference_dir, igblast_databases].join(',')
    validate_params(pathsToValidate)

    // process the sample sheet, then concat all fastq files in
    // the barcode directories & name based on the sample (not barcode)
    concatenated_files = parse_sample_sheet(fastq_dir, sample_sheet, params.barcode_dirs)
 
    // since we are now preparing different report, we need to group the files accordingly
    concatenated_files
        .groupTuple(by: 3)
        .set{nanocomp_input}   

    // then run nanocomp for QC
    nanocomp(nanocomp_input)

    // get the reference files
    Channel.fromPath("${reference_dir}/*.fasta")
        .map{ it -> [it.getSimpleName().replaceAll("_all_imgt_refs", ""), it]}
        .set{reference_files}

    // combine the reference files with the concatenated files tuple so that in the next
    // step we make sure we use the right reference file
    // adapted from https://stackoverflow.com/questions/72010030/how-combine-two-channels-by-tuple-elements-in-different-positions
    concatenated_files
        .map { tuple( it[2], *it ) }
        .combine ( reference_files, by: 0 )
        .map { it[1..-1] }
        .set {ab_selection_input}

    // select antibody reads (w/ minimap2 alignment)
    ab_reads = select_ab_reads(ab_selection_input)
    
    // trim the ends (by default polyA tail, user can specify if they have primers)
    trimmed_reads = read_trimming(ab_reads, trim_3p, trim_5p).trimmed_reads_renamed
    
    // annotation and grouping
    trimmed_read_fasta = fastq_to_fasta(trimmed_reads) //IgBLAST needs fasta not fastq
    annotation_grouping_pre_consensus(
        trimmed_read_fasta, 
        igblast_databases, 
        num_consensus)
    pre_consensus_table = annotation_grouping_pre_consensus.out.pre_consensus_table
    
    // get all of the read names + sequence ids into a channel
    annotation_grouping_pre_consensus.out.read_names_for_consensus
        .collect(flat: false) // ensures this step doesn't start until all files are ready
        .flatten() // gets rid of the extra list layer
        .map { it -> [it.baseName.minus("_read_names"), it] } // reformat to [sequence_id, file]
        .set{read_names_for_consensus}

    // do the same for the consensus starting points
    annotation_grouping_pre_consensus.out.starting_points_for_consensus
        .collect(flat: false) // ensures this step doesn't start until all files are ready
        .flatten() // gets rid of the extra list layer
        .map { it -> [it.baseName.minus("_starting_point_name"), it] } // reformat to [sequence_id, file]
        .set{starting_points_for_consensus}

    // join the read names and starting points for consensus, based on sequence id
    read_names_for_consensus
        .map { tuple( it[0], *it[1..-1] ) } // the asterisk unlists it 
        .combine(starting_points_for_consensus, by: 0)
        .set{read_names_and_starting_points_for_consensus}

    // now add the sample ID to this channel
    read_names_and_starting_points_for_consensus
        .map { it -> [it[0].replaceAll(~/_.._count_.*/, ""), *it] } // remove the chain, count etc leaving just sample ID
        .set {read_names_and_starting_points_for_consensus_with_sample_ID}

    // add the trimmed reads to the channel
    // first need to reshape the trimmed reads so that the sample ID is in the first position
    trimmed_reads
        .map { tuple( it[0][0], *it[1..-1] ) } // the asterisk unlists it 
        .set {trimmed_reads_with_sample_ID}

    // add the trimmed reads to the channel with read names and starting points
    read_names_and_starting_points_for_consensus_with_sample_ID
        .combine(trimmed_reads_with_sample_ID, by: 0)
        .set{input_for_consensus}

    // actually make the consensus 
    consensus_sequences = take_consensus(input_for_consensus, medaka_model)

    // annotation & grouping/analysis of consensus sequences
    // need to get the consensus sequences from [sequence_id, file] into [sample_id, [file1, file2, file3, etc]]
     consensus_sequences
        .map { it -> [it[0].replaceAll(~/_.._count_.*/, ""), *it[1..-1]] } // replace sequence ID with sample ID
        .groupTuple() // group by sample ID
        .set{consensus_sequences_by_sample_id}
     
    // make a channel of the metadata so that we can add it back in (we need it for the next step)
    trimmed_reads
        .map { tuple( it[0][0], it[0] ) } // this is just [sample ID, [metadata]]
        .set {sampleid_metadata}

    sampleid_metadata
        .join(consensus_sequences_by_sample_id)
        .map { it[1..-1] }
        .set{final_files_for_annotation}
    
    annotation_grouping_post_consensus(
        final_files_for_annotation,
        igblast_databases)
    
    // for now only using productive annotation in the report
    // i think it might confuse end users to give the unproductive sequences too?
    productive_only_annotation = annotation_grouping_post_consensus.out.only_productive_sequences

    // need to reformat this to combine with the nanocomp stuff
    // basically, get the metadata out of this structure: [[meta1, meta2, meta3], path]
    // into this structure: [meta3, meta2, meta1, path]
    // then can group by the report number
    productive_only_annotation
        .map{ it -> [it[0][2]] + [it[0][1]] + [it[0][0]] + it}
        .map{it -> it[0, 1, 2, 4]}
        .groupTuple(by: 0)
        .set{annotation_ordered_for_report}

    // for the nanocomp output, we need the report group first
    // can achieve this by just deleting the first two elements (sample name and species)
    // because they are already in the annotation_ordered_for_report channel and we will
    // join those two together
    nanocomp.out.txt
        .map{it -> it[2..-1]}
        .set{nanocomp_text_ordered_for_report}

    nanocomp.out.png
        .map{it -> it[2..-1]}
        .set{nanocomp_pics_ordered_for_report}

    // now do the joining 
    annotation_ordered_for_report
        .join(nanocomp_text_ordered_for_report)
        .join(nanocomp_pics_ordered_for_report)
        .set{reportable_data}
    
    versions()

    // make report
    //report_template = Channel.fromPath(params.report_template)
    //css_file = Channel.fromPath(params.css_file)
    //logo = Channel.fromPath(params.logo_file)
    
     report(
        version,
        sample_sheet,
        report_title,
        report_template,
        reportable_data,
        css_file,
        logo
    )
}