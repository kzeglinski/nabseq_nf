// adapted from the nf-core module: https://github.com/nf-core/modules/tree/master/modules/nf-core/minimap2/align
// (replaced meta$id with just meta) because sample name is our only metadata
// also set --secondary=no since we don't want secondary alignments and removed the bam output option (& samtools) since it was not needed 
// finally, changed to also emit the reads it aligned (because we need them in the next step of subsetting)
// idk why 

process minimap2_alignment {
    tag "$meta"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::minimap2=2.24' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/minimap2%3A2.24--h7132678_1' :
        'quay.io/biocontainers/minimap2:2.24--h7132678_1' }"

    input:
    tuple val(meta), path(reads)
    path reference

    output:
    tuple val(meta), path("*.paf"), path(reads), optional: true, emit: paf_reads
    tuple val(meta), path(reads), path(reference), path("*.paf"), optional: true, emit: for_racon
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta}"

    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        --secondary=no \\
        -x map-ont \\
        "${reference ?: reads}" \\
        "$reads" \\
        -o ${prefix}.paf


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
