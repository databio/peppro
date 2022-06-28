# Detailed installation instructions

This guide walks you through the minutiae of how to install each prerequisite component.  We'll presume you're installing this in a Linux environment.  If not the case, you'll need to go to each tool's respective site to find alternative installation approaches and options.

You have several options for installing the software prerequisites: 1) use a container, either [a single container](run-container.md) or with a [multi-container environment manager](run-bulker.md), in which case you need only either `docker` or `singularity`; 2) [install via `conda`](run-conda.md) or 3) install all prerequisites natively. We'll install everything natively in this guide. 

## 1. Install required software

**Python packages**. The pipeline uses [`pypiper`](http://pypiper.readthedocs.io/en/latest/) to run a single sample, [`looper`](http://looper.readthedocs.io/en/latest/) to handle multi-sample projects (for either local or cluster computation), [`pararead`](https://github.com/databio/pararead) for parallel processing sequence reads, [`refgenie`](http://refgenie.databio.org/en/latest/), optionally, to organize, build, and manage reference assemblies, [`cutadapt`](https://cutadapt.readthedocs.io/) to remove adapters, and a handful of `python` libraries. You can do a user-specific install of required `python` packages using the included `requirements.txt` file in the pipeline directory:  
```console
pip install --user -r requirements.txt
```

**Required executables**. We will need some common bioinformatics tools installed. The complete list (including optional tools) is specified in the pipeline configuration file (pipelines/peppro.yaml) tools section.
The following tools are used by the pipeline by default:  

* [bedtools (v2.30.0+)](http://bedtools.readthedocs.io/en/latest/)
* [bowtie2 (v2.4.2+)](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/) *installed via `pip`*
* [fastq-pair](https://github.com/linsalrob/fastq-pair.git)
* [flash](https://ccb.jhu.edu/software/FLASH/)
* [picard](https://broadinstitute.github.io/picard/) *which is required by the `python` package `pypiper`*
* [preseq (v2.0.3)](http://smithlabresearch.org/software/preseq/)
* [samtools (v1.14+)](http://www.htslib.org/)
* [seqkit](https://bioinf.shenwei.me/seqkit/)
* [seqtk](https://github.com/lh3/seqtk)
* Two UCSC tools (v3.5.1+)
    * [bigWigCat (v4)](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)
    * [wigToBigWig (v4)](https://www.encodeproject.org/software/wigtobigwig/)

We'll install each of these pieces of software before moving forward. Let's create an initial working directory to download and make all of this software.
```console
mkdir tools && cd tools/
```

### bedtools
We're going to install from source, but if you would prefer to install from a package manager, you can follow the instructions in the [bedtools' installation guide](https://bedtools.readthedocs.io/en/latest/content/installation.html).
```console
wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
tar -zxvf bedtools-2.30.0.tar.gz
rm bedtools-2.30.0.tar.gz
cd bedtools2
make
```

Now, let's add `bedtools` to our `PATH` environment variable.  Look here to [learn more about the concept of environment variables](https://www.digitalocean.com/community/tutorials/how-to-read-and-set-environmental-and-shell-variables-on-a-linux-vps) if you are unfamiliar.
```console
export PATH="$PATH:/path/to/tools/bedtools2/bin/"
```

### bowtie2
Next, let's install `bowtie2`.  For more more specific instruction, [read the author's installation guide](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#building-from-source).
```console
cd ../
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.4.2/bowtie2-2.4.2-source.zip
unzip bowtie2-2.4.2-source.zip
rm bowtie2-2.4.2-source.zip
cd bowtie2-2.4.2/
make
```
Again, let's add `bowtie2` to our `PATH` environment variable:
```
export PATH="$PATH:/path/to/tools/bowtie2-2.3.4.1/"
```
Great! On to the next one. 

### fastq_pair
Finally, because PRO-seq treats read1 differently than read2 in paired-end data, we need to re-sync paired-end files after processing.  We [use `fastq_pair`](https://github.com/linsalrob/fastq-pair/blob/master/INSTALLATION.md) to do so efficiently.
```console
cd ../
git clone https://github.com/linsalrob/fastq-pair.git
cd fastq-pair/
mkdir build
cd build/
cmake3 ..
make
make install
```

### flash

To obtain a plot to evaluate library quality when we have paired-end reads, [we use `FLASH`](https://ccb.jhu.edu/software/FLASH/) to generate a distribution of reads.  
```console
cd ../../
wget http://ccb.jhu.edu/software/FLASH/FLASH-1.2.11-Linux-x86_64.tar.gz
tar xvfz FLASH-1.2.11-Linux-x86_64.tar.gz
```

And let's add `FLASH` to our `PATH` environment variable:
```
export PATH="$PATH:/path/to/tools/FLASH-1.2.11-Linux-x86_64/"
```

### picard
`PEPPRO` is built using `PyPiper` and relies upon the `PyPiper NGSTK` tool kit which itself employs `Picard`.  [Read the `picard` installation guide](http://broadinstitute.github.io/picard/) for more assistance.
```console
wget https://github.com/broadinstitute/picard/releases/download/2.26.10/picard.jar
chmod +x picard.jar
```
Create an environmental variable pointing to the `picard.jar` file called `PICARD`.  Alternatively, [update the `peppro.yaml` file](https://github.com/databio/peppro/blob/master/pipelines/peppro.yaml) with the full `PATH` to the `picard.jar` file.
```
export PICARD="/path/to/tools/picard.jar"
```

### HTSlib
`PEPPRO` uses `samtools`, and `samtools` requires `HTSlib` internally. So first we'll install `HTSlib`.
```console
wget https://github.com/samtools/htslib/releases/download/1.14/htslib-1.14.tar.bz2
tar xvfj htslib-1.14.tar.bz2
cd htslib-1.14/
./configure
```

Alternatively, if you do not have the ability to install `HTSlib` to the default location, you can specify using the `--prefix=/install/destination/dir/` option.  [Learn more about the `--prefix` option here](http://samtools.github.io/bcftools/howtos/install.html). Otherwise, we will install to the default location.
```console
make
make install
```

### samtools
Next up, `samtools`.
```console
wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
tar xvfj samtools-1.14.tar.bz2
rm samtools-1.14.tar.bz2
cd samtools-1.14/
./configure
```
Alternatively, if you do not have the ability to install `samtools` to the default location, you can specify using the `--prefix=/install/destination/dir/` option.  [Learn more about the `--prefix` option here](http://samtools.github.io/bcftools/howtos/install.html). Otherwise, we will install to the default location.
```console
make
make install
```

### preseq
The pipeline uses `preseq` to calculate library complexity. Check out the author's [page for more instruction](https://github.com/smithlabcode/preseq). 
```console
wget https://github.com/smithlabcode/preseq/releases/download/v2.0.3/preseq_v2.0.3.tar.bz2
tar xvfj preseq_v2.0.3.tar.bz2
cd preseq/
make all SAMTOOLS_DIR=/path/to/tools/samtools-1.14
```
Add to `PATH`!
```console
export PATH="$PATH:/path/to/tools/preseq/"
```

### seqkit
Let's grab `seqkit` now.  Check out [the author's installation guide](https://github.com/shenwei356/seqkit#installation) for more instruction if necessary. 
```console
cd ../
wget https://github.com/shenwei356/seqkit/releases/download/v2.1.0/seqkit_linux_amd64.tar.gz
tar -zxvf seqkit_linux_amd64.tar.gz
```
And then make sure that executable is in our `PATH`.
```console
export PATH="$PATH:/path/to/tools/"
```

### UCSC utilities
Finally, we need a few of the UCSC utilities.  You can install the [entire set of tools](https://github.com/ENCODE-DCC/kentUtils) should you choose, but here we'll just grab the subset that we need.
```
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod 755 wigToBigWig
chmod 755 bigWigCat
```

That should do it!  

## 2. Install R packages

`PEPPRO` uses `R` to generate quality control plots.  These are technically **optional** and the pipeline will run without them, but you would not get any QC plots.  If you need to but don't have [R installed, you can follow these instructions](https://cran.r-project.org/doc/manuals/r-release/R-admin.html).  We'll use and install the necessary packages in this example.  Here is the list of required packages:

 - [data.table (v1.14.2)](https://cran.r-project.org/package=data.table)
 - [devtools (v2.4.3)](https://cran.r-project.org/web/packages/devtools/index.html)
 - [GenomicDistributions (v1.3.2)](http://code.databio.org/GenomicDistributions/index.html)
 - [ggplot2 (v3.3.5)](https://cran.r-project.org/package=ggplot2)
 - [pepr (v0.4.0)](http://code.databio.org/pepr/)
 - [optigrab (v0.9.2.1)](https://cran.r-project.org/web/packages/optigrab/index.html)

To install the needed packages, enter the following command in the pipeline folder:
```console
Rscript -e 'install.packages("devtools")'
Rscript -e 'devtools::install_github("pepkit/pepr")'
Rscript -e 'install.packages("BiocManager")'
Rscript -e 'BiocManager::install("GenomicRanges")'
Rscript -e 'devtools::install_github("databio/GenomicDistributions")'
wget "http://big.databio.org/GenomicDistributionsData/GenomicDistributionsData_0.0.2.tar.gz"
Rscript -e 'install.packages("GenomicDistributionsData_0.0.2.tar.gz", type="source", repos=NULL)'
Rscript -e 'devtools::install(file.path("PEPPROr/"), dependencies=TRUE, repos="https://cloud.r-project.org/")'
```

To extract files quicker, `PEPPRO` can also utilize `pigz` in place of `gzip` if you have it installed.  Let's go ahead and do that now. It's not required, but it can help speed everything up when you have many samples to process.
```
wget https://zlib.net/pigz/pigz-2.7.tar.gz
tar xvfz pigz-2.7.tar.gz
rm pigz-2.7.tar.gz
cd pigz-2.7/
make
cd ../
```
Don't forget to add this to your `PATH` too!
```
export PATH="$PATH:/path/to/tools/pigz-2.7/"
```

## 3. Download genomic assets using `refgenie`

`PEPPRO` can use [`refgenie`](http://refgenie.databio.org/) to simplify asset management for alignment, quality control reports, and some outputs. You can initialize a `refgenie` config file like this:

```console
pip install refgenie
export REFGENIE=genome_config.yaml
refgenie init -c $REFGENIE
```

Add the `export REFGENIE=genome_config.yaml` line to your `.bashrc` or `.profile` to ensure it persists. 

Next, pull the assets you need. Replace `hg38` in the example below if you need to use a different genome assembly. If these assets are not available automatically for your genome of interest, then you'll need to [build them](annotation.md). Download these required assets with this command:

```console
refgenie pull hg38/fasta hg38/bowtie2_index hg38/refgene_anno hg38/ensembl_gtf hg38/ensembl_rb
refgenie build hg38/feat_annotation
```
PEPPRO also requires `bowtie2_index` for any pre-alignment genomes:

```console
refgenie pull human_rDNA/fasta human_rDNA/bowtie2_index
```

That's it! Everything we need to run `PEPPRO` to its full potential should be installed.

## 4. Confirm installation 

You can confirm the pipeline is now executable natively using the included `checkinstall` script.  This can either be run directly from the `peppro/` repository...

```console
./checkinstall
```

or from the web:
```console
curl -sSL https://raw.githubusercontent.com/databio/peppro/checkinstall | bash
```
