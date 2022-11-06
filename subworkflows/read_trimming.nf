// trims reads and removes the rc from the read names

include { cutadapt } from '../modules/local/cutadapt'

process rename_cutadapt_output {
    tag "$sample_id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::seqkit=2.3.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0' }"
    
    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path('*_ab_reads_trimmed_renamed.fastq'), emit: reads

    script:
    """
    # extract the reads
    seqkit replace -p " rc" -r "" $reads > "${sample_id}_ab_reads_trimmed_renamed.fastq"

    """
}

workflow read_trimming {
    take:
        ab_reads
        trim_3p
        trim_5p
        
    main:
        // run cutadapt
        trimmed_reads = cutadapt(ab_reads, trim_3p, trim_5p).reads  
        
        // rename the trimmed reads
        trimmed_reads_renamed = rename_cutadapt_output(trimmed_reads).reads
    
    emit:
        trimmed_reads_renamed
}
