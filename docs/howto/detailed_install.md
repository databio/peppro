# Detailed installation instructions

This guide walks you through the minutiae of how to install each prerequisite component.

## 1: Clone the `PEPPRO` pipeline

To begin, we need to get the `PEPPRO` pipeline itself.  The pipeline is hosted on [github](https://github.com/databio/peppro).  If you don't have git installed, follow the [git installation instructions](https://git-scm.com/download/linux), and here is a [brief introduction to git](https://guides.github.com/introduction/git-handbook/). To install `PEPPRO`, you can use one of the following methods:

* using SSH:
```
git clone git@github.com:databio/peppro.git
```
* using HTTPS:
```
git clone https://github.com/databio/peppro.git
```

We'll use HTTPS in this example.  From an open terminal, let's first create a directory we'll use to run through this guide:
```
mkdir peppro_tutorial
```

Let's move into our newly created directory and create a few more folders that we'll use later.
```
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
Success! If you had any issues, feel free to [reach out to us with questions](https://github.com/databio/peppro/issues).  Otherwise, let's move on to installing additional software.

## 2: Install required software

You have two options for installing the software prerequisites: 1) use a container, in which case you need only either `docker` or `singularity`; or 2) install all prerequisites natively. 

To use `PEPPRO`, we need the following software:
**Python packages**. The pipeline uses [`pypiper`](http://pypiper.readthedocs.io/en/latest/) to run a single sample, [`looper`](http://looper.readthedocs.io/en/latest/) to handle multi-sample projects (for either local or cluster computation), [`pararead`](https://github.com/databio/pararead) for parallel processing sequence reads, [`refgenie`](http://refgenie.databio.org/en/latest/) to organize and build reference assemblies, and the common `python` libraries [`numpy`](https://www.numpy.org/) and [`pandas`](https://pandas.pydata.org/). You can do a user-specific install using the included requirements.txt file in the pipeline directory:  
```
pip install --user -r requirements.txt
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
* UCSC tools (v3.5.1)
    * [wigToBigWig (v4)](https://www.encodeproject.org/software/wigtobigwig/)
    * [bigWigCat (v4)](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/)

We'll install each of these pieces of software before moving forward.  Let's start right at the beginning and install `bedtools`.  We're going to install from source, but if you would prefer to install from a package manager, you can follow the instructions in the [bedtools' installation guide](http://bedtools.readthedocs.io/en/latest/content/installation.html).
```
cd tools/
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz
rm bedtools-2.25.0.tar.gz
cd bedtools2
make
```

Now, let's add `bedtools` to our `PATH` environment variable.  Look here to [learn more about the concept of environment variables](https://www.digitalocean.com/community/tutorials/how-to-read-and-set-environmental-and-shell-variables-on-a-linux-vps) if you are unfamiliar.
```
export PATH="$PATH:/path/to/peppro_tutorial/tools/bedtools2/bin/"
```
Next, let's install `bowtie2`.
```
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

* [TODO] seqkit
* [TODO] fastp
* [TODO] preseq

Next up, `samtools`.
```
wget https://sourceforge.net/projects/samtools../files/samtools/1.9/samtools-1.9.tar.bz2
tar xvfj samtools-1.9.tar.bz2
rm samtools-1.9.tar.bz2
cd samtools-1.9
/configure
```
Alternatively, if you do not have the ability to install `samtools` to the default location, you can specify using the `--prefix=/install/destination/dir/` option.  [Learn more about the `--prefix` option here](http://samtools.github.io/bcftools/howtos/install.html).
```
make
make install
```
As for our other tools, add `samtools` to our `PATH` environment variable:
```
export PATH="$PATH:/path/to/peppro_tutorial/tools/samtools-1.9/"
```

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
That should do it!  Now we'll [install some **optional** packages](../tutorial.md#install-optional-software).  Of course, these are not required, but for the purposes of this tutorial we're going to be completionists.

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

Before we analyze anything, we also need a reference genome.  `PEPPRO` uses `refgenie` genomes.  For the purposes of this tutorial, we'll just download pre-built genomes.  Follow the `'refgenie` instructions if you'd like to [build your own reference genome](https://github.com/databio/refgenie). First, let's change into our `genomes/` folder.
```
cd /path/to/peppro_tutorial/genomes/
wget http://big.databio.org/refgenomes/hg38.tgz
wget http://cloud.databio.org.s3.amazonaws.com/refgenomes/human_repeats_170502.tgz
wget http://cloud.databio.org.s3.amazonaws.com/refgenomes/rCRSd_170502.tgz
tar xvfz hg38.tgz
tar xvfz human_repeats_170502.tgz
tar xvfz rCRSd_170502.tgz
rm hg38.tgz
rm human_repeats_170502.tgz
rm rCRSd_170502.tgz
```

## 4: Point the pipeline to your Refgenie assemblies

Let's also create another environment variable that points to our genomes.
```
export GENOMES="/path/to/peppro_tutorial/genomes/
```
(Don't forget to add this to your `.bashrc` or `.profile` to ensure it persists).


## 5: Download or create annotation files

To calculate TSS enrichments, you will need a [TSS annotation file](http://big.databio.org/refgenomes/) in your reference genome directory.  If a pre-built version for your genome of interest isn't present, you can quickly create that file yourself. In the reference genome directory, you can perform the following commands for in this example, `hg38`:
```
wget -O hg38_TSS_full.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz \
zcat hg38_TSS_full.txt.gz | \
  awk  '{if($4=="+"){print $3"\t"$5"\t"$5"\t"$4"\t"$13}else{print $3"\t"$6"\t"$6"\t"$4"\t"$13}}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_TSS.tsv
```

We also have [downloadable pre-built genome annotation files](http://big.databio.org/peppro/) for `hg38`, `hg19`, `mm10`, and `mm9` that you can use to annotate the reads and peaks.  These files annotate 3' and 5' UTR, Exonic, Intronic, Intergenic, Promoter, and Promoter Flanking Regions of the corresponding genome as indicated in Ensembl or UCSC.  Simply move the corresponding genome annotation file into the `peppro/anno` folder.  Once present in the `peppro/anno` folder you don't need to do anything else as the pipeline will look there automatically.   Alternatively, you can use the `--anno-name` pipeline option to directly point to this file when running.  You can also [learn how to create a custom annotation file](annotation_files.md) to calculate coverage using your own features of interest.

Alright! Time to setup the pipeline configuration files and run our sample.
