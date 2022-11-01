// convert fastq to fasta (needed for IgBLAST)
process fastq_to_fasta {
    tag {sample_name}
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::seqkit=2.3.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit%3A2.3.1--h9ee0642_0' :
        'quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0' }"

    input:
    tuple val(sample_name), path(reads)

    output:
    tuple val(sample_name), path("${sample_name}_ab_reads_trimmed.fasta")

    script:
    """
    seqkit fq2fa --threads $task.cpus $reads -o "${sample_name}_ab_reads_trimmed.fasta"
    """

}
