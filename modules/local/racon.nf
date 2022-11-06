// adapted from the nf-core racon module https://github.com/nf-core/modules/blob/master/modules/nf-core/racon/main.nf

process racon {
    tag "$sequence_id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::racon=1.4.20" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/racon:1.4.20--h9a82719_1' :
        'quay.io/biocontainers/racon:1.4.20--h9a82719_1' }"

    input:
    tuple val(sequence_id), path(reads), path(assembly), path(paf)

    output:
    tuple val(sequence_id), path(reads), path('*_racon_consensus.fasta') , emit: results

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sequence_id}"
    """
    racon -t "$task.cpus" \\
        "${reads}" \\
        "${paf}" \\
        $args \\
        -w 5000 -u -g -8 -x -6 -m 8 \\
        "${assembly}" > \\
        ${prefix}_racon_consensus.fasta

    """
}