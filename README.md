![NAb-seq](./modules/report/nabseq_logo.png)
# NAb-seq (nextflow version)
NAb-seq is an accurate, rapid and cost-effective method for antibody long-read sequencing in hybridoma 
cell lines and single B cells. 

# Useful links
* The [NAb-seq manuscript](https://www.tandfonline.com/doi/full/10.1080/19420862.2022.2106621)
* The [NAb-seq documentation website](https://kzeglinski.github.io/nab-seq/index.html)
* The [original version of NAb-seq](https://github.com/kzeglinski/nabseq_old)
* The [data from the NAb-seq manuscript](https://www.ebi.ac.uk/ena/browser/view/PRJEB51442?show=reads) 

# Usage
Also see the [documentation](https://kzeglinski.github.io/nab-seq/index.html) for more usage information
```
Usage: nextflow run ./nabseq_nf/main.nf --fastq_dir [input path] --organism [organism name] --sample_sheet [sample sheet]
--help                : prints this help message

Required arguments:
--out_dir             : where the output files will be written to (default: "$projectDir/results)
--fastq_dir           : where the input fastq files are located
--organism            : name of the organism whose reference sequences to use (default: "rat")
--sample_sheet        : location of the .csv sample sheet (format: barcode01,sample_x)

Optional (only needed for advanced users)
--num_consensus       : maximum number of consensus sequences to generate for each sample (default: 999)
--igblast_databases   : location of the igblast databases (default: "$projectDir/references/igblast/")
--reference_sequences : location of the reference sequences (default: "$projectDir/references/reference_sequences/")
--trim_3p             : pattern for cutadapt to trim off the 3' end (default: "A{20}N{90}")
--trim_5p             : pattern for cutadapt to trim off the 5' end (default: "N{90}T{20}")
--medaka_model        : model to use for medaka (depends on your basecalling model, default:"r941_min_sup_g507")
--report_title        : title to use for naming the report (default: "NAb-seq report")
```

## Input
NAb-seq requires as input:
* The path of a directory where it should write the results
* The path of the directory containing the basecalled and demultiplexed nanopore reads
* The name of the organism whose reference sequences should be used (rat and mouse are built-in, and you can create custom references for any organism you like 
* A sample sheet (CSV) in the format: barcode01,sample_x

## Output
NABseq will produce a report containing sample/run information, QC metrics and the annotation results of productive sequences from all samples. All consensus sequences in FASTA format, as well as intermediate files (NanoComp QC results, IgBLAST results) are written to the results directory also


