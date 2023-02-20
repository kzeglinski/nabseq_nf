// make the report

process write_report {
    tag "creating report ${report_number}"
    label 'process_low'
    publishDir "${params.out_dir}/report", mode: 'copy', pattern: "*.html"
    stageInMode 'copy'

    conda (params.enable_conda ? 'r::r-tidyverse=1.2.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse%3A1.2.1' :
        'quay.io/biocontainers/r-tidyverse:1.2.1' }"

    input:
    val version
    path sample_sheet
    val report_title
    path report_template
    tuple val(report_number), val(organism), val(sample_name), path(productive_only_annotation), path(nanocomp_txt), path(nanocomp_pngs)
    path css_file
    path logo

    output:
    path('*.html'), emit: report

    script:

    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(rmarkdown)
    rmarkdown::render(
        "$report_template",
        params = list(
            version = "$version",
            sample_sheet = "$sample_sheet",
            report_title = "${report_title}_${report_number}", 
            annotation_table = "$productive_only_annotation",
            nanocomp_pngs = "$nanocomp_pngs",
            nanocomp_txt = "$nanocomp_txt",
            report_number = $report_number,
            css_file = "$css_file"),
            output_file = "${report_title}_${report_number}.html"
        )

    """
}

workflow report  {
    take:
        version
        sample_sheet
        report_title
        report_template
        reportable_data
        css_file
        logo
        
    main:
        // process this in R
        write_report(version,
        sample_sheet,
        report_title,
        report_template,
        reportable_data,
        css_file,
        logo)
        report = write_report.out.report

    emit:
        report
        

}