// this workflow trims the 3p and 5p ends of reads (by default, polyA and polyT but you can specify your own)
// adapted from the nf-core module https://github.com/nf-core/modules/tree/master/modules/nf-core/cutadapt

// have updated the trimming strategy a bit to account for (1) forgetting to trim adapters/barcodes
// and (2) those weird cases where polyT and polyA appear in the middle of a read. might be a little
// aggressive but i guess thats better than the consensus sequences being cooked 


process cutadapt {
    tag "$meta"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::cutadapt=3.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cutadapt:3.4--py39h38f01e4_1' :
        'quay.io/biocontainers/cutadapt:3.4--py39h38f01e4_1' }"

    input:
    tuple val(meta), path(reads)
    val trim_3p
    val trim_5p

    output:
    tuple val(meta), path('*_ab_reads_trimmed.fastq'), emit: reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    cutadapt \\
        --cores $task.cpus \\
        $args \\
        -a $trim_3p \\
        -g $trim_5p \\
        -n 4 \\
        -m 300 \\
        -o "${meta}_ab_reads_trimmed.fastq" \\
        $reads
    """
}