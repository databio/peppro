# Detailed installation instructions

This guide walks you through the minutiae of how to install each prerequisite component.  We'll presume you're installing this in a Linux environment.  If not the case, you'll need to go to each tool's respective site to find alternative installation approaches and options.

## 1: Clone the `PEPPRO` pipeline

To begin, we need to get the `PEPPRO` pipeline itself.  The pipeline is hosted on [github](https://github.com/databio/peppro).  If you don't have git installed, follow the [git installation instructions](https://git-scm.com/download/linux), and here is a [brief introduction to git](https://guides.github.com/introduction/git-handbook/). 

From an open terminal, let's first create a directory we'll use to run through this guide:
```console
mkdir peppro_tutorial
```

Let's move into our newly created directory and create a few more folders that we'll use later.
```console
cd peppro_tutorial/
mkdir data
mkdir genomes
mkdir processed
mkdir templates
mkdir tools
cd tools/
```

Time to get PEPPRO!
```
git clone https://github.com/databio/peppro.git
```
(You could instead use SSH instead of HTTPS with `git clone git@github.com:databio/peppro.git`). Success! If you had any issues, feel free to [reach out to us with questions](https://github.com/databio/peppro/issues).  Otherwise, let's move on to installing additional software.

## 2: Install required software

You have two options for installing the software prerequisites: 1) use a container, in which case you need only either `docker` or `singularity`; or 2) install all prerequisites natively. We'll install everything natively in this guide. If you want to try the container approach, read [PEPPRO in containers](container.md).

To use `PEPPRO`, we need the following software:
**Python packages**. The pipeline uses [`pypiper`](http://pypiper.readthedocs.io/en/latest/) to run a single sample, [`looper`](http://looper.readthedocs.io/en/latest/) to handle multi-sample projects (for either local or cluster computation), [`pararead`](https://github.com/databio/pararead) for parallel processing sequence reads, [`refgenie`](http://refgenie.databio.org/en/latest/) to organize and build reference assemblies, [`cutadapt`](https://cutadapt.readthedocs.io/) to remove adapters for single-end data or optionally for paired-end, and the common `python` libraries [`numpy`](https://www.numpy.org/) and [`pandas`](https://pandas.pydata.org/). You can do a user-specific install using the included requirements.txt file in the pipeline directory:  
```console
pip install --user -r requirements.txt
```
Remember to add your user specific install location to your `PATH`.
```console
export PATH="$PATH:$HOME/.local/bin/"
```

**Required executables**. We will need some common bioinformatics tools installed. The complete list (including optional tools) is specified in the pipeline configuration file (pipelines/peppro.yaml) tools section.
The following tools are used by the pipeline:  

* [bedtools (v2.25.0+)](http://bedtools.readthedocs.io/en/latest/)
* [bowtie2 (v2.2.9+)](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
* [samtools (v1.7)](http://www.htslib.org/)
* [seqkit](https://bioinf.shenwei.me/seqkit/)
* [fastp](https://github.com/OpenGene/fastp)
* [seqtk](https://github.com/lh3/seqtk)
* [preseq](http://smithlabresearch.org/software/preseq/)
* [fastq-pair](https://github.com/linsalrob/fastq-pair.git)
* [picard](http://broadinstitute.github.io/picard/)
* UCSC tools (v3.5.1)
    * [wigToBigWig (v4)](https://www.encodeproject.org/software/wigtobigwig/)
    * [bigWigCat (v4)](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)

#### bedtools
We'll install each of these pieces of software before moving forward.  Let's start right at the beginning and install `bedtools`.  We're going to install from source, but if you would prefer to install from a package manager, you can follow the instructions in the [bedtools' installation guide](http://bedtools.readthedocs.io/en/latest/content/installation.html).
```console
cd tools/
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
rm bedtools-2.25.0.tar.gz
cd bedtools2
make
```

Now, let's add `bedtools` to our `PATH` environment variable.  Look here to [learn more about the concept of environment variables](https://www.digitalocean.com/community/tutorials/how-to-read-and-set-environmental-and-shell-variables-on-a-linux-vps) if you are unfamiliar.
```console
export PATH="$PATH:/path/to/peppro_tutorial/tools/bedtools2/bin/"
```

#### bowtie2
Next, let's install `bowtie2`.  For more more specific instruction, [read the author's installation guide](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#building-from-source).
```console
cd ../
wget https://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.3.4.1/bowtie2-2.3.4.1-source.zip
unzip bowtie2-2.3.4.1-source.zip
rm bowtie2-2.3.4.1-source.zip
cd bowtie2-2.3.4.1
make
```
Again, let's add `bowtie2` to our `PATH` environment variable:
```
export PATH="$PATH:/path/to/peppro_tutorial/tools/bowtie2-2.3.4.1/"
```
Great! On to the next one. 

#### samtools
Next up, `samtools`.
```console
wget https://sourceforge.net/projects/samtools../files/samtools/1.9/samtools-1.9.tar.bz2
tar xvfj samtools-1.9.tar.bz2
rm samtools-1.9.tar.bz2
cd samtools-1.9
./configure
```
Alternatively, if you do not have the ability to install `samtools` to the default location, you can specify using the `--prefix=/install/destination/dir/` option.  [Learn more about the `--prefix` option here](http://samtools.github.io/bcftools/howtos/install.html).
```console
make
make install
```
As for our other tools, add `samtools` to our `PATH` environment variable:
```
export PATH="$PATH:/path/to/peppro_tutorial/tools/samtools-1.9/"
```

#### seqkit
Let's grab `seqkit` now.  Check out [the author's installation guide](https://github.com/shenwei356/seqkit#installation) for more instruction if necessary. 
```console
cd ../
wget https://github.com/shenwei356/seqkit/releases/download/v0.10.1/seqkit_linux_amd64.tar.gz
tar -zxvf seqkit_linux_amd64.tar.gz
```
And then make sure that executable is in our `PATH`.
```console
export PATH="$PATH:/path/to/peppro_tutorial/tools/"
```

#### fastp
Next on our list is `fastp`. Check out their [install instructions](https://github.com/OpenGene/fastp#get-fastp) if necessary.
```console
git clone https://github.com/OpenGene/fastp.git
cd fastp/
make
make install
```
Add to `PATH`!
```console
export PATH="$PATH:/path/to/peppro_tutorial/tools/fastp/"
```

#### preseq
The pipeline uses `preseq` to calculate library complexity. Check out the author's [page for more instruction](https://github.com/smithlabcode/preseq).
```console
cd ../
wget http://smithlabresearch.org/downloads/preseq_linux_v2.0.tar.bz2
tar xvfj preseq_linux_v2.0.tar.bz2
```
Add to `PATH`!
```console
export PATH="$PATH:/path/to/peppro_tutorial/tools/preseq_v2.0/"
```

#### fastq_pair
Finally, because PRO-seq treats read1 differently than read2 in paired-end data, we need to resync paired-end files after processing.  We [use `fastq_pair`](https://github.com/linsalrob/fastq-pair/blob/master/INSTALLATION.md) to do so efficiently.
```console
git clone https://github.com/linsalrob/fastq-pair.git
cd fastq-pair/
mkdir build
cd build/
cmake3 ..
make
make install
```

#### picard
`PEPPRO` is built using `PyPiper` and relies upon the `PyPiper NGSTK` tool kit which itself employs `Picard`.  [Read the `picard` installation guide](http://broadinstitute.github.io/picard/) for more assistance.
```console
cd ../../
wget https://github.com/broadinstitute/picard/releases/download/2.20.3/picard.jar
chmod +x picard.jar
```
Create an environmental variable pointing to the `picard.jar` file called `$PICARD`.  Alternatively, [update the `peppro.yaml` file](https://github.com/databio/peppro/blob/master/pipelines/peppro.yaml) with the full PATH to the `picard.jar` file.
```
export PICARD="/path/to/peppro_tutorial/tools/picard.jar"
```

#### UCSC utilities
Finally, we need a few of the UCSC utilities.  You can install the [entire set of tools](https://github.com/ENCODE-DCC/kentUtils) should you choose, but here we'll just grab the subset that we need.
```
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigCat
chmod 755 wigToBigWig
chmod 755 bigWigCat
```
Add our `tools/` directory to our `PATH` environment variable.
```
export PATH="$PATH:/path/to/peppro_tutorial/tools/"
```
That should do it!  Now we'll install some **optional** packages.  Of course, these are not required, but for the purposes of this tutorial we're going to be completionists.

### Optional software

`PEPPRO` uses `R` to generate quality control plots.  These are **optional** and the pipeline will run without them, but you would not get any QC plots.  If you need to don't have [R installed, you can follow these instructions](https://cran.r-project.org/doc/manuals/r-release/R-admin.html).  We'll use and install the necessary packages in this example.  Here is the list of required packages:

 - [data.table (v1.11.2)](https://cran.r-project.org/package=data.table)
 - [GenomicDistributions (v0.5)](http://code.databio.org/GenomicDistributions/index.html)
 - [ggplot2 (v2.2.1)](https://cran.r-project.org/package=ggplot2)
 - [pepr (v0.2.1)](http://code.databio.org/pepr/)
 - [optigrab (v0.9.2.1)](https://cran.r-project.org/web/packages/optigrab/index.html)

To install the needed packages, enter the following command in the pipeline folder:
```
Rscript -e 'install.packages("PEPPROr", repos=NULL, type="source")'
```
To extract files quicker, `PEPPRO` can also utilize `pigz` in place of `gzip` if you have it installed.  Let's go ahead and do that now. It's not required, but it can help speed everything up when you have many samples to process.
```
cd /path/to/peppro_tutorial/tools/
wget http://zlib.net/pigz/pigz-2.4.tar.gz
tar xvfz pigz-2.4.tar.gz
rm pigz-2.4.tar.gz
cd pigz-2.4/
make
```
Don't forget to add this to your `PATH` too!
```
export PATH="$PATH:/path/to/peppro_tutorial/tools/pigz-2.4/"
```
That's it! Everything we need to run `PEPPRO` to its full potential should be installed.  If you are interested and have experience using containers, you can check out the [alternate installation methods](../install.md#121-use-containers).

### Create environment variables

We also need to create some environment variables to help point `looper` to where we keep our data files and our tools.  You may either set the environment variables up, like we're going to do now, or you may simply hard code the necessary locations in our configuration files.
First, let's create a `PROCESSED` variable that represents the location where we want to save output.
```
export PROCESSED="/path/to/peppro_tutorial/processed/"
```
Second, we'll create a variable representing the root path to all our tools named `CODEBASE`.
```
export CODEBASE="/path/to/peppro_tutorial/tools/"
```
(Add these environment variables to your `.bashrc` or `.profile` so you don't have to always do this step).
Fantastic! Now that we have the pipeline and its requirements installed, we're ready to get our reference genome(s).


## 3: Download a reference genome


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

PEPPRO also uses [refgenie](https://refgenie.databio.org) to manage a variety of annotation files for quality control plots. Downloading them is very easy:

```
refgenie pull -g hg38 -a ensembl_gtf tss_annotation feat_annotation pi_tss pi_body
```
Replace `hg38` if you need to use a different genome assembly. If these assets are not available automatically for your genome of interest, then you'll need to [build them](annotation.md).


Alright! Time to setup the pipeline configuration files and run our sample.
