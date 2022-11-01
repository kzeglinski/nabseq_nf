//adapted from the nf-core NanoPlot module https://nf-co.re/modules/nanoplot
process nanocomp {
    tag "nanocomp_on_all_reads"
    label 'process_low'
    publishDir "${params.out_dir}/nanocomp", mode: 'copy'

    conda (params.enable_conda ? 'bioconda::nanocomp=1.19.3' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanocomp%3A1.19.3--pyhdfd78af_0' :
        'quay.io/biocontainers/nanocomp:1.19.3--pyhdfd78af_0' }"

    input:
    path(ontfile)

    output:
    path("*.html"), emit: html
    path("*.png") , emit: png
    path("*.txt") , emit: txt
    path("*.log") , emit: log
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    NanoComp \\
        $args \\
        -t $task.cpus \\
        --fastq $ontfile
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanocomp: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
