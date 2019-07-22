# Getting started

## 1: Clone the `PEPPRO` pipeline

Clone the pipeline:
```
git clone https://github.com/databio/peppro.git
```

## 2: Install required software

`PEPPRO` requires a series of publicly-available, common bioinformatics tools including: [samtools](http://www.htslib.org/), [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [seqkit](https://bioinf.shenwei.me/seqkit/), [fastp](https://github.com/OpenGene/fastp), [seqtk](https://github.com/lh3/seqtk), [preseq](http://smithlabresearch.org/software/preseq/), [fastq-pair](https://github.com/linsalrob/fastq-pair.git), [picard](http://broadinstitute.github.io/picard/), [wigToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/), and [bigWigCat](http://hgdownload.soe.ucsc.edu/admin/exe/).

### Python packages

`PEPPRO` uses several packages under the hood. Make sure you're up-to-date with a user-specific install:

```{bash}
cd peppro
pip install --user -r requirements.txt
```

### R package

`PEPPRO` uses R to produce QC plots, and we include an R package for these functions. From the `peppro/` directory:
```console
Rscript -e 'install.packages("PEPPROr", repos=NULL, type="source")'
```

### Optional software

Optionally, `PEPPRO` can mix and match tools for adapter removal, read trimming, deduplication, and reverse complementation.  The use of `fqdedup`, in particular, is useful if you wish to minimize memory use at the expense of speed.  We suggest using the default tools simply due to the fact that `fastx toolkit` has not been supported since 2012.

`seqOutBias` can be used to take into account the mappability at a given read length to filter the sample signal.

*Optional tools:* [fqdedup](https://github.com/guertinlab/fqdedup), [fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/), [seqOutBias](https://github.com/guertinlab/seqOutBias), [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc), and [pigz (v2.3.4+)](https://zlib.net/pigz/).

## 3: Download `refgenie` assemblies

The pipeline relies on [`refgenie` assemblies](http://refgenie.databio.org/en/dev/install/) for alignment.  First, initialize a folder for genome indexes and the `refgenie` config file.

```console
export REFGENIE=your_genome_folder/genome_config.yaml
refgenie init -c $REFGENIE
```

Then, just pull the assets you need.

```console
refgenie pull -g hg38 -a bowtie2
refgenie pull -g rCRSd -a bowtie2
refgenie pull -g human_repeats -a bowtie2
```

(Add `REFGENIE` to your .bashrc or .profile to ensure it persists). Alternatively, you can skip the `REFGENIE` variable and simply change the value of the `resources.genome_config` option in the [`pipeline_config.yaml`](https://github.com/databio/peppro/blob/master/pipelines/peppro.yaml) file to point to the folder where you stored the assemblies. 

## 4: Run the pipeline script directly

The pipeline at its core is just a python script, and you can run it on the command line for a single sample (see [command-line usage](usage)), which you can also get on the command line by running `pipelines/peppro.py --help`. You just need to pass a few command-line parameters to specify sample name, reference genome, input files, etc. Here's the basic command to run the included small test example through the pipeline:

```console
/pipelines/peppro.py \
  --sample-name test \
  --genome hg38 \
  --input examples/data/test_r1.fq.gz \
  --single-or-paired single \
  -O $HOME/peppro_example/
```

This test example takes less than 5 minutes to complete. Read more about how to [run the test sample using `Looper`](howto/use_looper.md) with the included [example `peppro_test.yaml` file](https://github.com/databio/peppro/blob/master/examples/meta/peppro_test.yaml).

# 5. Next steps

This is just the beginning. For your next step, take a look at one of these user guides:

- [Extended tutorial for running a single sample](tutorial.md)
- [Running on multiple samples with looper](howto/use_looper.md)
- [Running the pipeline directly in a container](howto/use_container.md)
- See other detailed user guide links in the side menu
