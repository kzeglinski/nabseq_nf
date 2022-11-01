// the two-step consensus calling process (racon, then medaka)

include { racon } from '../modules/local/racon' 
include { medaka } from '../modules/local/medaka' 
include { minimap2_alignment } from '../modules/local/minimap2' 

process prepare_consensus_inputs{
    tag "$sequence_id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::seqkit=2.3.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0' }"
    
    input:
    tuple val(sample_id), val(sequence_id), path(read_names), path(starting_point_name), path(trimmed_reads)

    output:
    tuple val(sequence_id), path("${sequence_id}_reads.fastq"), emit: reads
    path("${sequence_id}_starting_point.fastq"), emit: starting_point

    script:
    """
    # extract the reads
    seqkit grep --threads ${task.cpus} -n -r -f $read_names $trimmed_reads | seqkit replace -p "\\_.*" -r "" | seqkit rename > "${sequence_id}_reads.fastq"

    # extract starting copy
    seqkit grep --threads ${task.cpus} -n -r -f $starting_point_name $trimmed_reads | seqkit replace -p ".*" -r $sequence_id > "${sequence_id}_starting_point.fastq"
    """
}

workflow take_consensus {
    take:
        consensus_input
        medaka_model
        
    main:
        // use the read name files to make fasta/fastq of the actual reads 
        prepare_consensus_inputs(consensus_input)
        consensus_input_reads = prepare_consensus_inputs.out.reads
        consensus_starting_point = prepare_consensus_inputs.out.starting_point
        
        // run racon
        minimap2_alignment(consensus_input_reads, consensus_starting_point)
        racon(minimap2_alignment.out.for_racon)

        // run medaka
        medaka(racon.out.results, medaka_model)
        consensus_sequences = medaka.out.consensus

    emit:
        consensus_sequences
}
