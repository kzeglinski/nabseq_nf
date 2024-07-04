// adapted from the nf-core medaka module https://github.com/nf-core/modules/blob/master/modules/nf-core/medaka/main.nf

process medaka {
    tag "$sequence_id"
    label 'process_high'
    publishDir "${params.out_dir}/consensus_sequences", mode: 'copy'


    conda (params.enable_conda ? "bioconda::medaka=1.11.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.11.3--py39h05d5c5e_0' :
        'quay.io/biocontainers/medaka:1.11.3--py39h05d5c5e_0' }"

    input:
    tuple val(sequence_id), path(reads), path(assembly)
    val (medaka_model)

    output:
    tuple val(sequence_id), path("*.fasta"), emit: consensus

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sequence_id}"
    """
    medaka_consensus \\
        -t $task.cpus \\
        $args \\
        -i $reads \\
        -d $assembly \\
        -m $medaka_model \\
        -o ./
    mv consensus.fasta ${prefix}.fasta

    """
}
