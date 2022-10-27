// module for running IgBLAST
// general layout is based on the nf-core modules
process igblast {
    tag "$meta"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::igblast=1.19.0' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/igblast%3A1.19.0--pl5321h3928612_0' :
        'quay.io/biocontainers/igblast:1.19.0--pl5321h3928612_0' }"
    
    input:
    tuple val(meta), path(reads)
    val organism
    val igblast_databases
    env IGDATA
    env IGBLASTDB

    output:
    tuple val(meta), path('*_pre_consensus_igblast.tsv'), emit: airr_table

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    # run igblast
    # outfmt 19 = AIRR format (tsv, easy to use in downstream steps)
    igblastn -germline_db_V ${igblast_databases}/databases/imgt_${organism}_V \
        -germline_db_J ${igblast_databases}/databases/imgt_${organism}_J \
        -germline_db_D ${igblast_databases}/databases/imgt_${organism}_D \
        -c_region_db ${igblast_databases}/databases/imgt_${organism}_C \
        -organism $organism \
        -query $reads \
        -auxiliary_data ${igblast_databases}/igdata/optional_file/${organism}_gl.aux \
        -show_translation \
        -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -num_alignments_C 1 \
        -outfmt 19 > ${meta}_pre_consensus_igblast.tsv
    """
}