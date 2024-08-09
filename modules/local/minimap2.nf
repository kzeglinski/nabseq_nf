// adapted from the nf-core module: https://github.com/nf-core/modules/tree/master/modules/nf-core/minimap2/align
// (replaced meta$id with just meta) because sample name is our only metadata
// update as of august 2024: now we are using sam output instead of paf (scales better)

process minimap2_alignment {
    tag "$prefix"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::minimap2=2.24' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2%3A2.24--h7132678_1' :
        'quay.io/biocontainers/minimap2:2.24--h7132678_1' }"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.sam"), path(reads), optional: true, emit: sam_reads
    tuple val(meta), path(reads), path(reference), path("*.sam"), optional: true, emit: for_racon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // allow for a bunch of metadata (although the first element should be sample name)
    if(meta instanceof Collection) {
        prefix = meta[0]
    } else {
        prefix = meta
    }

    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        --secondary=no \\
        --sam-hit-only \\
        -ax map-ont \\
        "${reference ?: reads}" \\
        "$reads" \\
        -o ${prefix}.sam

    """
}
