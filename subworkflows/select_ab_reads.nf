// do minimap2 alignment to the specified reference
// then, use seqkit to select the reads that align

include { minimap2_alignment } from '../modules/local/minimap2' 

process subset_aligned_reads {
    tag {prefix}
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::seqkit=2.3.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(paf_file), path(reads)

    output:
    tuple val(meta), path("*_ab_reads.fastq")

    script:

    // allow for a bunch of metadata (although the first element should be sample name)
    if(meta instanceof Collection) {
        prefix = meta[0]
    } else {
        prefix = meta
    }

    """
    # aligned read names are the first column of the PAF file
    cut -f1 $paf_file > "${prefix}_ab_read_names.txt"
    
    # use these names as a pattern for seqkit grep to find
    seqkit grep --by-name --use-regexp --threads ${task.cpus} \
    -f "${prefix}_ab_read_names.txt" \
    $reads \
    -o "${prefix}_ab_reads.fastq"
    """

}

workflow select_ab_reads {
    take:
        ab_selection_input
    
    main:
        // split input into reads + meta data and then reference file
        ab_selection_input
            .map{it[-1]}
            .set{reference}

        ab_selection_input
            .map{it[0..-2]}
            .map{tuple([it[0], it[2], it[3]], it[1])}
            .set{minimap2_input}

        // align reads using minimap2
        minimap2_alignment(minimap2_input, reference)

        // subset the aligned reads using seqkit
        ab_reads = subset_aligned_reads(minimap2_alignment.out.paf_reads)
        
    emit:
        ab_reads

}