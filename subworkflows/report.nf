// make the report

process write_report {
    tag "creating report ${report_number}"
    label 'process_low'
    publishDir "${params.out_dir}/report", mode: 'copy', pattern: "*.zip"
    stageInMode 'copy'
    container "library://kzeglinski/nabseq/nabseq-report:v0.0.3"

    input:
    val report_number
    path sample_sheet
    path report_template
    path versions
    tuple val(report_number), val(organism), val(sample_name), path(productive_only_annotation), path(nanocomp_htmls), path(flag_data)

    output:
    path('*.zip'), emit: report

    script:

    """
    export DENO_DIR="\$PWD"
    export XDG_CACHE_HOME="/tmp/quarto_cache_home"
    export XDG_DATA_HOME="/tmp/quarto_data_home"

    unzip template_files.zip
    quarto render report_template.qmd --to html --output "NAb-seq report ${report_number}.html"
    mkdir NAb-seq_report_${report_number}
    cp "NAb-seq report ${report_number}.html" NAb-seq_report_${report_number}
    cp favicon.svg NAb-seq_report_${report_number}
    cp nabseq_logo.png NAb-seq_report_${report_number}
    cp -r report_template_files/ NAb-seq_report_${report_number}
    zip -r NAb-seq_report_${report_number}.zip NAb-seq_report_${report_number}/
    """
}

workflow report  {
    take:
        sample_sheet
        reportable_data
        report_template
        versions

    main:
        reportable_data
            .map{it -> it[0]}
            .set{report_number}

        // process this in R
        write_report(
            report_number,
            sample_sheet,
            report_template,
            versions,
            reportable_data)
        report = write_report.out.report

    emit:
        report
}