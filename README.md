# PEPPRO

PEPPRO is a pipeline designed to process PRO-seq data. For more information see: http://code.databio.org/PEPPRO/

## Required software

PEPPRO uses a series of publicly-available, common bioinformatics tools including:

* [samtools](http://www.htslib.org/)
* [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [seqkit](https://bioinf.shenwei.me/seqkit/)
* [fastp](https://github.com/OpenGene/fastp)
* [seqtk](https://github.com/lh3/seqtk)
* [seqOutBias](https://github.com/guertinlab/seqOutBias)
  
## Optional software

Alternatively, `PEPPRO` can mix and match tools for adapter removal, read trimming, deduplication, and reverse complementation.  The use of `fqdedup`, in particular, is useful if you wish to minimize memory use at the expense of increased speed.  I also suggest using the default required tools simply due to the fact that the `fastx toolkit` has not been supported since 2012 and issues with reads utilizing newer Phred quality scores can cause problems.

*Optional tools:*
* [fqdedup](https://github.com/guertinlab/fqdedup)
* [cutadapt](https://cutadapt.readthedocs.io/)
* [fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)

## Contributing

Pull requests welcome. Active development should occur in a development or feature branch.

## Contributors

* Jason Smith, jasonsmith@virginia.edu
* Nathan Sheffield, nathan@code.databio.org
* Mike Guertin, mjg7y@virginia.edu
* Others... (add your name)
