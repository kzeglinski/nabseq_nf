// do minimap2 alignment to the specified reference
// then, use samtools to subset

include { minimap2_alignment } from '../modules/local/minimap2'

process subset_aligned_reads {
    tag {prefix}
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::samtools=1.20' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_1' :
    'quay.io/biocontainers/samtools:1.20--h50ea8bc_1' }"

    input:
    tuple val(meta), path(sam_file), path(reads)

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
    samtools fastq -@ ${task.cpus} $sam_file > ${prefix}_ab_reads.fastq
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
        ab_reads = subset_aligned_reads(minimap2_alignment.out.sam_reads)

    emit:
        ab_reads

}