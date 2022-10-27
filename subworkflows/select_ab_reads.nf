// do minimap2 alignment to the specified reference
// then, use seqkit to select the reads that align

include { minimap2_alignment} from '../modules/local/minimap2' 

process subset_aligned_reads {
    tag {sample_name}

    conda (params.enable_conda ? 'bioconda::seqkit=2.3.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0' }"

    input:
    tuple val(sample_name), path(paf_file), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}_ab_reads.fastq")

    script:
    """
    # aligned read names are the first column of the PAF file
    cut -f1 $paf_file > "${sample_name}_ab_read_names.txt"
    
    # use these names as a pattern for seqkit grep to find
    seqkit grep --by-name --use-regexp \
    -f "${sample_name}_ab_read_names.txt" \
    $reads \
    -o "${sample_name}_ab_reads.fastq"
    """

}

workflow select_ab_reads {
    take:
        concatenated_files
        reference_sequence
    
    main:
        // align reads using minimap2
        minimap2_alignment(concatenated_files, reference_sequence)

        // subset the aligned reads using seqkit
        ab_reads = subset_aligned_reads(minimap2_alignment.out.paf_reads)

    emit:
        ab_reads

}