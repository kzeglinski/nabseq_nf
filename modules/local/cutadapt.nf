// this workflow trims the 3p and 5p ends of reads (by default, polyA and polyT but you can specify your own)
// adapted from the nf-core module https://github.com/nf-core/modules/tree/master/modules/nf-core/cutadapt

// have updated the trimming strategy a bit to also re-orient the reads. This shouldn't have been a problem
// but sometimes it was creating some really cooked consensus sequences with constant regions at both ends
// re-orienting seems to have fixed this problem


process cutadapt {
    tag "$prefix"
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
    def trim_pattern = "-a ${trim_3p}"
    if (!trim_5p.isEmpty()) {
        trim_pattern = trim_pattern.concat(" -g ${trim_5p}")
    }

    // allow for a bunch of metadata (although the first element should be sample name)
    if(meta instanceof Collection) {
        prefix = meta[0]
    } else {
        prefix = meta
    }

    """
    cutadapt \\
        --cores $task.cpus \\
        $args \\
        $trim_pattern \\
        -n 2 \\
        -m 300 \\
        --revcomp \\
        -o "${prefix}_ab_reads_trimmed.fastq" \\
        $reads
    """
}