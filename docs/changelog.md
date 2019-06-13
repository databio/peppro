# Change log
All notable changes to this project will be documented in this file.

## [0.7.2] -- 2019-06-13

### Changed
- Fix --new-start bug
- Fix container libgsl library
- Adjust seqOutBias usage and flexibility

## [0.7.1] -- 2019-06-12

### Changed
- Update to new refgenie usage
- Update container

## [0.7.0] -- 2019-06-11

### Changed
- Fix multiple input file handling
- Update preseq calculation and plotting and handle small samples
- Fix TSS plotting and score calculation
- Include bedtools in required software

### Added
- Add project level library complexity summarizer
- Add docs
- Add container
- Add small test example

## [0.6.0] -- 2019-06-04

### Changed
- Use PEPPRO R package for QC plotting and analysis
- Default uses intermediate read files to produce a deduplicated and non-deduplicated aligned BAM file
- Check for cutadapt version for multicore processing
- Update bigWig production to be variable step formatted wig

### Added
- Handle GRO-seq data as input
- Created and included a PEPPRO R package for standard functions
- Add preseq requirement and plot library complexity curve
- Calculate and report read depth
- Calculate and report the percent of adapter contamination
- Produce fraction of reads in features plot

## [0.5.1] -- 2019-05-09

### Changed
- Simplify and clarify prealignment steps

### Added
- Perform a pre-check that all required tools are callable
