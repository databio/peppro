# Install and run <img src="../img/peppro_logo.svg" alt="PEPPRO" class="img-fluid" style="max-height:50px; margin-top:-15px; margin-bottom:-10px">

## 1: Clone the `PEPPRO` pipeline

```
git clone https://github.com/databio/peppro.git
```

## 2: Install required software

PEPPRO requires a set of Python and R packages to run.

### Python packages

`PEPPRO` uses several packages under the hood. Make sure you're up-to-date with a user-specific install:

```{bash}
cd peppro
pip install --user -r requirements.txt
```

### R package

`PEPPRO` uses R to produce QC plots, and we include an R package for these functions.  The `PEPPRO` package relies on a handful of additional packages. 

To install the prerequisite packages from the command line:
```console
Rscript -e 'install.packages("devtools")'
Rscript -e 'devtools::install_github("pepkit/pepr")'
Rscript -e 'install.packages("BiocManager")'
Rscript -e 'BiocManager::install("GenomicRanges")'
Rscript -e 'devtools::install_github("databio/GenomicDistributions")'
```

Then, install the `PEPPRO` package. From the `peppro/` directory:
```console
Rscript -e 'devtools::install(file.path("PEPPROr/"), dependencies=TRUE, repos="https://cloud.r-project.org/")'

```

### Tools

The pipeline also relies on a set of publicly available bioinformatic tools, but if you don't want to install the prerequisite software used by PEPPRO natively, you can learn to [run PEPPRO using containers](container.md) and skip this step.

Otherwise, you'll need to install the following: [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html), [bigWigCat](http://hgdownload.soe.ucsc.edu/admin/exe/), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [fastq-pair](https://github.com/linsalrob/fastq-pair.git), [flash](https://ccb.jhu.edu/software/FLASH/), [picard](https://broadinstitute.github.io/picard/), [preseq](http://smithlabresearch.org/software/preseq/), [seqkit](https://bioinf.shenwei.me/seqkit/), [samtools](http://www.htslib.org/), [seqtk](https://github.com/lh3/seqtk), and [wigToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/). If you need help, we have [detailed installation instructions](detailed_install.md) for installing these.

## 3: Download `refgenie` assets

PEPPRO uses [`refgenie`](http://refgenie.databio.org/) assets for alignment. If you haven't already, initialize a refgenie config file like this:

```console
export REFGENIE=your_genome_folder/genome_config.yaml
refgenie init -c $REFGENIE
```

Add the `export REFGENIE` line to your `.bashrc` or `.profile` to ensure it persists. 

Next, pull the assets you need. Replace `hg38` in the example below if you need to use a different genome assembly. If these assets are not available automatically for your genome of interest, then you'll need to [build them](annotation.md). Download these required assets with this command:

```console
refgenie pull -g hg38 -a bowtie2_index ensembl_gtf ensembl_rb refgene_anno feat_annotation 
```
PEPPRO also requires `bowtie2_index` for any pre-alignment genomes:

```console
refgenie pull -g human_rDNA -a bowtie2_index
```

### Optional software

Optionally, `PEPPRO` can mix and match tools for adapter removal, read trimming, deduplication, and reverse complementation.  The use of `fqdedup`, in particular, is useful if you wish to minimize memory use at the expense of speed.  We suggest using the default tools simply due to the fact that `fastx toolkit` has not been supported since 2012. `seqOutBias` can be used to take into account the mappability at a given read length to filter the sample signal.

*Optional tools:*

* [fastp](https://github.com/OpenGene/fastp)
* [fqdedup](https://github.com/guertinlab/fqdedup)
* [fastx toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
* [seqOutBias](https://github.com/guertinlab/seqOutBias)
* [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc)
* [pigz (v2.3.4+)](https://zlib.net/pigz/)

## 4: Run an example project through PEPPRO

Start by running the example project (`peppro_test.yaml`) in the [`examples/meta/`](https://github.com/databio/peppro/tree/master/examples/meta) folder. PEPPRO uses a project management tool called [looper](https://looper.databio.org) to run the pipeline across samples in a project. Let's use the `-d` argument to do a *dry run*, which will create job scripts for every sample in a project, but will not execute them:

```
cd peppro
looper run -d examples/meta/peppro_test.yaml
```

If the looper executable is not in your `$PATH`, add the following line to your `.bashrc` or `.profile`:
```
export PATH=$PATH:~/.local/bin
```
If that worked, let's actually run the example by taking out the `-d` flag:

```console
looper run examples/meta/peppro_test.yaml
```

Or, if you're using [`bulker`](https://bulker.databio.org/en/latest/) to run the pipeline in containers:

```console
bulker activate databio/peppro
looper run examples/meta/peppro_test.yaml
```

There are lots of other cool things you can do with looper, like dry runs, summarize results, check on pipeline run status, clean intermediate files to save disk space, lump multiple samples into one job, and more. For details, consult the [`looper` docs](http://looper.databio.org/).

## 5: Configure your project files

To run your own samples, you'll need to organize them in **PEP format**, which is explained in [how to create a PEP](https://pepkit.github.io/docs/home/) and is universal to all pipelines that read PEPs, including `PEPPRO`. To get you started, there are examples you can adapt in the `examples/` folder (*e.g.* [example test PEP](https://github.com/databio/peppro/tree/master/examples/meta/peppro_test.yaml)). In short, you need two files for your project:

  1. project config file -- describes output locations, pointers to data, etc.
  2. sample annotation file -- comma-separated value (CSV) list of your samples.

The sample annotation file must specify these columns:

- sample_name
- library (*e.g.* 'PRO', 'PROSEQ', 'PRO-seq', 'GRO', 'GROSEQ', 'GRO-seq')
- organism (*e.g.* 'human' or 'mouse')
- read1
- read2 (if paired)
- anything else you wish to include

## Next steps

This is just the beginning. For your next step, the [extended tutorial](tutorial.md) will walk you through a real project. Or, take a look at one of other detailed user guide links in the side menu.
