# Changelog
All notable changes to this project will be documented in this file.

## [v0.2.0]

First release of the nextflow version of NAb-seq. Key changes from the original version are:
- it's now written in nextflow (duh)
- References for mouse and rat now included in NAb-seq by default
    Previously, users had to manually create their own reference sequences and
    IgBLAST databases. Now, we provide these for mouse and rat, so you don't
    need to do it yourself. If you would like to create your own custom reference,
    please see the reference sequence guide (in the navbar at the top of the
    page)
- Annotation of constant regions is now performed with IgBLAST
    The original version of NAb-seq was written before IgBLAST could be used on
    constant regions, so constant annotation was performed using minimap2. To
    simplify the bioinformatics, we have now moved to using IgBLAST.
    **Please note that this doesn't affect the results in any way, and we provide pre-built IgBLAST C gene databases for mouse and rat.**
- Improved trimming & reorientation of reads
    We have now implemented an improved trimming strategy which also reorients
    the reads. This produces cleaner consensus sequences, with less 'junk'
    (adapters, polyA tails) hanging off the end. There is also now support for
    custom adapter trimming, by changing the default values of the `trim_3p`
    and `trim_5p` parameters (see below for more information).

## [v0.2.1]

Small fix: corrected the logic for choosing a starting copy for consensus calling. Now, we choose the longest read with a complete V(D)J region (previously it was just the longest read, so in rare cases the consensus sequence generated would be incomplete)

## [v0.2.2]

New features:
- you can now specify a `species` column in the sample sheet (so that hybridomas/b cells of different species can be processed in the same run)
- you can now group samples into reports using a `group` column in the sample sheet (useful for pooling many samples from different labs and issuing a separate report to each one)

## [v0.2.3]

Bug fixes:
- added `.fq` and `.fq.gz` to the list of valid fastq file extensions
- corrected a problem where, when running NAb-seq with many samples, the pipeline would crash when reshaping channels for consensus calling
- corrected a bug where consensus starting copies would not be written for groups without a complete VDJ, leading to the consensus calling process silently failing
- updated the README to show correct sample sheet structure

New features:
- added ability to run NAb-seq on a single directory full of fastqs rather than having to make separate folders for each barcode

## [v0.2.4]

- added a new subworkflow to collect the versions of all the tools in modules into a versions.yml file which will be stored in the output directory.
- added validation of the input paths

## [v0.3.0]
Exciting update!!!!!
- reports are now interactive!
- you can download all of the tables, and the sequences in FASTA format
- also introduced a flagging system to classify hybridomas/cells as 'normal', containing 'multiple H and or L chains' or 'missing H and or L chain'
- added the --trim_to_antibody option (on by default) which trims the sequences to only contain from the start of the V region to the end of the C region
- changed the default value of the --trim_5p parameter to be the sequence of the PCB111 TSO, seeing as the PCB109 kit is now discontinued
- new 'versions' workflow creates a file with all software versions. these are incorportated into the reports as well, in a short paragraph explaining the NAb-seq method
