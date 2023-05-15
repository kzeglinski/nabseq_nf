// ADAPTED FROM https://stackoverflow.com/questions/74039553/nextflow-rename-barcodes-and-concatenate-reads-within-barcodes

process concat_reads {

    tag { sample_name }
    label 'process_low'
    publishDir "${params.out_dir}/concat_reads", mode: 'copy', failOnError: true

    input:
    tuple val(sample_name), path(fastq_files), val(species), val(report_group)

    output:
    tuple val(sample_name), path("${sample_name}.${extn}"), val(species), val(report_group)

    script:
    if( fastq_files.every { it.name.endsWith('.fastq.gz') } )
        extn = 'fastq.gz'
    else if( fastq_files.every { it.name.endsWith('.fastq') } )
        extn = 'fastq'
    else if( fastq_files.every { it.name.endsWith('.fq.gz') } )
        extn = 'fq.gz'
    else if( fastq_files.every { it.name.endsWith('.fq') } )
        extn = 'fq'
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
        barcode_dirs
    
    main:
        fastq_extns = [ '.fastq', '.fastq.gz' , '.fq', '.fq.gz' ] // update to cover all possible fastq file extensions
        
        if( barcode_dirs ){ // deal with the case that fastqs are located in folders named by barcode
            Channel.fromPath(sample_sheet)
                | splitCsv(header: true)
                | map{row ->
                    full_path = fastq_dir + "/" + "${row.barcode}"
                    all_files = file(full_path).listFiles()
                    fastq_files = all_files.findAll { fn ->
                        fastq_extns.find { fn.name.endsWith( it ) }
                    }
                    tuple(row.sample_name, fastq_files, row.species, row.report_group)
                }
                | concat_reads
                | set{concatenated_file_tuple} 
        }
        else { // deal with the case that fastqs are all located in the same folder and named barcode01{something}.fq.gz
            Channel.fromPath(sample_sheet)
                | splitCsv(header: true)
                | map{row ->
                    fastq_files = file(fastq_dir).listFiles().findAll { fn ->
                        fn.name.startsWith( "${row.barcode}" ) && fastq_extns.find { fn.name.endsWith( it ) }
                    }
                    tuple(row.sample_name, fastq_files, row.species, row.report_group)
                }
                | concat_reads
                | set{concatenated_file_tuple} 
        }

    emit:
        concatenated_file_tuple

}