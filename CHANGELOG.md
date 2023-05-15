# Changelog
All notable changes to this project will be documented in this file.

## [v0.2.0]

First release of the nextflow version of NAb-seq

## [v0.2.1]

Small fix: corrected the logic for choosing a starting copy for consensus calling. Now, we choose the longest read with a complete V(D)J region (previously it was just the longest read, so in rare cases the consensus sequence generated would be incomplete)

## [v0.2.2]

New features: 
- you can now specify a `species` column in the sample sheet (so that hybridomas/b cells of different species can be processed in the same run)
- you can now group samples into reports using a `group` column in the sample sheet (useful for pooling many samples from different labs and issuing a separate report to each one)

## [v0.2.3]
Bug fixes:
- added `.fq` and `.fq.gz` to the list of valid fastq file extensions 

New features:
- 
