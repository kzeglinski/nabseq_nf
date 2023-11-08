//adapted from the nf-core NanoPlot module https://nf-co.re/modules/nanoplot
process nanocomp {
    tag "$sample_name"
    label 'process_low'
    publishDir "${params.out_dir}/nanocomp", mode: 'copy'

    conda (params.enable_conda ? 'bioconda::nanocomp=1.19.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocomp%3A1.19.3--pyhdfd78af_0' :
        'quay.io/biocontainers/nanocomp:1.19.3--pyhdfd78af_0' }"

    input:
    tuple val(sample_name), file(fastq_files), val(species), val(report_group)

    output:
    tuple val(sample_name), val(species), val(report_group), path("*.html"), emit: html

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    NanoComp \\
        $args \\
        -t $task.cpus \\
        --fastq $fastq_files \\
        --prefix "report_${report_group}_"

    """
}
