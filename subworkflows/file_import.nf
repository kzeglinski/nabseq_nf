// ADAPTED FROM https://stackoverflow.com/questions/74039553/nextflow-rename-barcodes-and-concatenate-reads-within-barcodes

process concat_reads {

    tag { sample_name }

    publishDir "${params.out_dir}/concat_reads", mode: 'copy'

    input:
    tuple val(sample_name), path(fastq_files)

    output:
    tuple val(sample_name), path("${sample_name}.${extn}")

    script:
    if( fastq_files.every { it.name.endsWith('.fastq.gz') } )
        extn = 'fastq.gz'
    else if( fastq_files.every { it.name.endsWith('.fastq') } )
        extn = 'fastq'
    else
        error "Concatentation of mixed filetypes is unsupported"

    """
    cat ${fastq_files} > "${sample_name}.${extn}"
    """
}

workflow parse_sample_sheet {
    take:
        fastq_dir
        sample_sheet
    
    main:
        fastq_extns = [ '.fastq', '.fastq.gz' ]
        Channel.fromPath(sample_sheet)
            |splitCsv()
            | map { dir, sample_name ->
                full_path = fastq_dir + "/" + dir // prepend the directory our barcoded directories live in
                all_files = file(full_path).listFiles()
                fastq_files = all_files.findAll { fn ->
                    fastq_extns.find { fn.name.endsWith( it ) }
                }
                tuple( sample_name, fastq_files )
            }
            | concat_reads
            | set{concatenated_file_tuple}
    
    emit:
        concatenated_file_tuple

}