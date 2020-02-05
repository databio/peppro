# Change log
All notable changes to this project will be documented in this file.

## [0.8.9] -- 2020-02-05

### Changed
- Corrected Pct_reads_too_short calculation

## [0.8.8] -- 2020-02-04

### Changed
- Updated degradation ratio calculation for PE data
- Corrected Pct_reads_too_short to be percent not fraction

## [0.8.7] -- 2020-02-04

### Changed
- Updated degradation ratio calculation for SE data
- Fixed report_fastq to properly handle SE, non UMI data

## [0.8.6] -- 2020-01-28

### Changed
- Update FRiF calculation to optionally follow a priority ranked method
- Update how adapter insertion distributions are plotted to be the same for SE or PE data
- Make cutadapt the default for adapter removal
- Streamline the use of Refgenie assets
    - Refgenie manages pause indicies
    - Refgenie manages feature annotations
    - Refgenie manages assets for mRNA contamination
    - Refgenie manages seqOutBias required suffixerator indicies
- Change pause index and mRNA contamination plots to histograms

### Added
- Add PRiF plot
- Require FLASH tool
- Produce sample level gene counts file as output
- Generate project level counts table including all samples X gene counts
- Report degradation metric for library quality
- Add BiocProject integration

## [0.8.1] -- 2019-07-15

### Changed
- Fix fraction in feature calculation
- Fix library complexity calculation for PE data
- Require fastq_pair tool
- Require cutadapt for SE processing

### Added
- Add QC plot for adapter insertion distribution

## [0.8.0] -- 2019-07-10

### Changed
- Update pause index calculation and required annotation files
- Update mRNA contamination required annotation files
- Move fastq processing to separate function
- Change handling of PE data
- Update TSS profiling

### Added
- Add modified fastq_pair tool to handle PE data properly
- Add pause index plotting
- Add mRNA contamination calculation and plotting
- Add fragment length distribution plotting for PE data

## [0.7.3] -- 2019-06-13

### Changed
- Fix missing gt in container
- Fix mappability bug with new refgenie

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
