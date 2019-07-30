# Install and run PEPPRO

## 1: Clone the `PEPPRO` pipeline

```
git clone https://github.com/databio/peppro.git
```
## 2: Download `refgenie` assets

PEPPRO uses [`refgenie`](http://refgenie.databio.org/) assets for alignment. If you haven't already, initialize a refgenie config file like this:

```console
pip install --user refgenie
export REFGENIE=your_genome_folder/genome_config.yaml
refgenie init -c $REFGENIE
```

Add the `export REFGENIE` line to your `.bashrc` or `.profile` to ensure it persists. Then, pull the assets you need. By default, that's these for human:

```console
refgenie pull -g hg38 -a bowtie2_index
refgenie pull -g human_rDNA -a bowtie2_index
refgenie pull -g rCRSd -a bowtie2_index
```

## 3: Install required software

If you don't want to install the prerequisite software used by PEPPRO, you can follow our tutorial on [running PEPPRO directly in a container](howto/use_container.md) and then skip this step. If you want to run it natively, you'll need to install the following: [samtools](http://www.htslib.org/), [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html), [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), [seqkit](https://bioinf.shenwei.me/seqkit/), [fastp](https://github.com/OpenGene/fastp), [seqtk](https://github.com/lh3/seqtk), [preseq](http://smithlabresearch.org/software/preseq/), [fastq-pair](https://github.com/linsalrob/fastq-pair.git), [wigToBigWig](http://hgdownload.soe.ucsc.edu/admin/exe/), and [bigWigCat](http://hgdownload.soe.ucsc.edu/admin/exe/).


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

*Optional tools:*

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

Or, if you're using containers, adjust the `--compute` argument accordingly:

```console
looper run examples/meta/peppro_test.yaml --compute docker
looper run examples/meta/peppro_test.yaml --compute singularity
```

There are lots of other cool things you can do with looper, like dry runs, summarize results, check on pipeline run status, clean intermediate files to save disk space, lump multiple samples into one job, and more. For details, consult the [`looper` docs](http://looper.databio.org/).

## 5: Configure your project files

To run your own samples, you'll need to organize them in **PEP format**, which is explained in [how to create a PEP](https://pepkit.github.io/docs/home/) and is universal to all pipelines that read PEPs, including `PEPPRO`. To get you started, there are examples you can adapt in the `examples/` folder (*e.g.* [example test PEP](https://github.com/databio/peppro/tree/master/examples/meta/peppro_test.yaml)). In short, you need two files for your project:

  1. project config file -- describes output locations, pointers to data, etc.
  2. sample annotation file -- comma-separated value (CSV) list of your samples.

The sample annotation file must specify these columns:

- sample_name
- library ('PRO' or 'PROSEQ' or 'PRO-seq')
- organism (e.g. 'human' or 'mouse')
- read1
- read2 (if paired)
- whatever else you want

## Next steps

This is just the beginning. For your next step, take a look at one of other detailed user guide links in the side menu.
