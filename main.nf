#!/usr/bin/env nextflow

// TODO: intro/copyright notice

// TODO: help message

// TODO: parameter validation
fastq_dir = params.fastq_dir
sample_sheet = params.sample_sheet

// TODO: introductory message that prints out parameters



/*
 * Bring in modules
 */
include { parse_sample_sheet } from './modules/local/file_import'
include { nanocomp } from './modules/local/nanocomp'

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

    // trim the 3' and 5' ends (by default polyA tail, user can specify if they have primers)

    // run igblast on the trimmed reads

    // group reads based on igblast output 

    // consensus calling 

    // re-igblast of consensus sequences

    // process this new output

    // make report

    // finished! print some kind of message to user

    // also email if they want?
}