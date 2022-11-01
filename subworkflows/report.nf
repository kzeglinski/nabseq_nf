// make the report

process write_report {
    tag "reporting_all_samples"
    label 'process_low'
    publishDir "${params.out_dir}/report", mode: 'copy', pattern: "*.html"
    stageInMode 'copy'

    conda (params.enable_conda ? 'r::r-tidyverse=1.2.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-tidyverse%3A1.2.1' :
        'quay.io/biocontainers/r-tidyverse:1.2.1' }"

    input:
    val version
    val organism
    path sample_sheet
    val report_title
    path report_template
    path productive_only_annotation
    path full_annotation
    path nanocomp_output
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
            version = "$version", organism = "$organism", 
            sample_sheet = "$sample_sheet",
            report_title = "$report_title", 
            annotation_table = "$productive_only_annotation",
            nanocomp_output = "$nanocomp_output",
            css_file = "$css_file"),
            output_file = "${report_title}.html"
        )

    """
}

workflow report  {
    take:
        version
        organism
        sample_sheet
        report_title
        report_template
        productive_only_annotation
        full_annotation
        nanocomp_output
        css_file
        logo
        
    main:
        // process this in R
        write_report(version,
        organism,
        sample_sheet,
        report_title,
        report_template,
        productive_only_annotation,
        full_annotation,
        nanocomp_output,
        css_file,
        logo)
        report = write_report.out.report

    emit:
        report
        

}