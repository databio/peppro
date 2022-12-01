# <img src="../img/peppro_logo.svg" alt="PEPPRO" class="img-fluid" style="max-height:50px; margin-top:-15px; margin-bottom:-10px"> pipeline step-by-step guide

In this guide, we'll walk you through the step by step procedure of running a tutorial PRO-seq dataset through the pipeline.  The output from this process is the same as you see in the [example PRO-seq output](browse_output.md) we've provided.  To use this tutorial, you should have a basic familiarity with [working in a command line driven environment](http://matt.might.net/articles/basic-unix/). You also need to have already installed `PEPPRO` prerequisites, which you can do following one of the [installation instructions](install.md).

# Tutorial using `refgenie` managed genome assets

## 1. Set up folders 

From an open terminal, let's first create a directory we'll use to run through this guide:
```console
mkdir peppro_tutorial
```

Let's point an environment variable to our tutorial location (change to match your local path) to link our tutorial samples with your local environment.
```console
export TUTORIAL=/path/to/your/peppro_tutorial
```

Let's move into our newly created directory and create a few more folders that we'll use later.
```console
cd peppro_tutorial/
mkdir data
mkdir processed
mkdir templates
mkdir tools
cd tools/
git clone https://github.com/databio/peppro.git
```

## 2: Initialize `refgenie` and download assets

`PEPPRO` can utilize [`refgenie`](http://refgenie.databio.org/) assets to simplify setup and sample analysis. Because assets are user-dependent, these files must always exist outside of any container system or alongside a native installation. Therefore, we still need to [install and initialize a refgenie config file.](http://refgenie.databio.org/en/latest/install/). For example:

```console
pip install refgenie
export REFGENIE=${TUTORIAL}/refgenie_config.yaml
refgenie init -c $REFGENIE
```

Add the `export REFGENIE` line to your `.bashrc` or `.profile` to ensure it persists. 

Next, pull the assets you need. Replace `hg38` in the example below if you need to use a different genome assembly. If these assets are not available automatically for your genome of interest, then you'll need to [build them](annotation.md).

```console
refgenie pull hg38/fasta hg38/bowtie2_index hg38/refgene_anno hg38/ensembl_gtf hg38/ensembl_rb
refgenie build hg38/feat_annotation
```

`PEPPRO` also requires a `fasta` and `bowtie2_index` asset for any pre-alignment genomes:

```console
refgenie pull human_rDNA/fasta human_rDNA/bowtie2_index
```

## 3. Download tutorial read files

We're going to work with some files a little larger than the test data included in the pipeline so we can see all the features included in a full run of the pipeline.  Go ahead and download the [tutorial_r1.fastq.gz](http://big.databio.org/peppro/fastq/tutorial_r1.fq.gz) and [tutorial_r2.fq.gz](http://big.databio.org/peppro/fastq/tutorial_r2.fq.gz) files. 
```console
wget http://big.databio.org/peppro/fastq/tutorial_r1.fq.gz
wget http://big.databio.org/peppro/fastq/tutorial_r2.fq.gz
```

To simplify the rest of this tutorial, let's put those files in a standard location we'll use for the rest of this guide. 
```console
mv tutorial_r1.fq.gz peppro/examples/data/
mv tutorial_r2.fq.gz peppro/examples/data/
```

## 4. Configure project files

We're going to use `looper` to analyze our data.  For that, we need to pass looper a configuration file.  This project config file describes your project. See [`looper` docs](https://looper.readthedocs.io/en/latest/) for details. A configuration file has been provided for you in the pipeline itself, conveniently named [`tutorial_refgenie.yaml`](https://github.com/databio/peppro/blob/master/examples/meta/tutorial_refgenie.yaml).  This configuration file also points to our sample.  In this case, we've provided a sample for you with the pipeline.  You don't have to do anything else at this point and may [skip right to running the sample if you'd like](tutorial.md#3-using-looper-to-run-the-pipeline).  Otherwise, we'll briefly touch on what those configuration files look like.

You can open the configuration file in your favorite text editor if you'd like to look closer.  For the purposes of the tutorial you may safely move past this step should you choose.
```console
cd peppro/examples/meta/
nano tutorial_refgenie.yaml
```
The following is what you should see in that configuration file.
```console
# Run tutorial samples through PEPPRO
name: PEPPRO_tutorial

pep_version: 2.0.0
sample_table: tutorial.csv

looper:
  output_dir: "${TUTORIAL}/processed/peppro/tutorial" 
  pipeline_interfaces: "${TUTORIAL}/tools/peppro/project_pipeline_interface.yaml"]

sample_modifiers:
  append:
    pipeline_interfaces: "${TUTORIAL}/tools/peppro/sample_pipeline_interface.yaml"
  derive:
    attributes: [read1, read2]
    sources:
        R1: "${TUTORIAL}/tools/peppro/examples/data/{sample_name}_r1.fq.gz"
        R2: "${TUTORIAL}/tools/peppro/examples/data/{sample_name}_r2.fq.gz"
  imply:
    - if:
        organism: ["human", "Homo sapiens", "Human", "Homo_sapiens"]
      then:
        genome: "hg38"
        prealignment_names: ["human_rDNA"]
```
There is also a sample annotation file referenced in our configuration file.  The sample annotation file contains metadata and other information about our sample. Just like before, this file, named [`tutorial.csv`](https://github.com/databio/peppro/blob/master/examples/meta/tutorial.csv) has been provided.  You may check it out if you wish, otherwise we're all set.

If you choose to open `tutorial.csv`, you should see the following:
```console
sample_name,organism,protocol,read_type,read1,read2
tutorial,human,PROSEQ,paired,R1,R2
```
That's it! Let's analyze that sample!

## 5. Use `looper` to run the pipeline
Looper requires a few variables and configuration files to work for the specific user. Let's get those set up now. `Looper` uses [`divvy`](https://divvy.databio.org/) to manage computing resource configuration so that projects and pipelines can easily travel among environments. For more detailed information, [check out the `looper` docs](https://looper.readthedocs.io/en/latest/running-on-a-cluster/). Let's set it up.
```console
cd $TUTORIAL
touch compute_config.yaml
```
Open that file in your favorite text editor.  We'll add in the following example for running locally.  You'll need to edit this file further for your own setup and you can [learn more about that in the `looper` docs](https://looper.readthedocs.io/en/latest/index.html).
```console
nano compute_config.yaml
```
Paste the following into compute_config.yaml
```console
adapters:
  CODE: looper.command
  JOBNAME: looper.job_name
  CORES: compute.cores
  LOGFILE: looper.log_file
  TIME: compute.time
  MEM: compute.mem

compute_packages:
  default:
    submission_template: templates/localhost_template.sub
    submission_command: .

```
Now, let's close and save that file and create an environment variable pointing to our configuration file.
```console
export DIVCFG="${TUTORIAL}/compute_config.yaml"
```

(Remember to add `DIVCFG` to your `.bashrc` or `.profile` to ensure it persists).
The `looper` environment configuration file points to submission template(s) in order to know how to run a samples locally or using cluster resources.  If you'd like to learn more, check out the [`DIVCFG` configuration file and submission templates](https://divvy.databio.org/). We're going to simply setup a local template for the purposes of this tutorial.  You can also easily create [templates for cluster or container use as well](https://github.com/pepkit/divcfg/tree/master/templates)!
Let's change to our `templates/` directory to make our first submission template.
```console
cd ${TUTORIAL}/templates/
nano localhost_template.sub
```             
Paste the following into the localhost_template.sub:
```console
#!/bin/bash

echo 'Compute node:' `hostname`
echo 'Start time:' `date +'%Y-%m-%d %T'`

{
{CODE}
} | tee {LOGFILE} --ignore-interrupts
```

Save and close that file, and return to the pipeline repository directory.
```console
cd ${TUTORIAL}/tools/peppro/
```

Now, we'll use `looper` to run the sample pipeline locally.
```console
looper run examples/meta/tutorial_refgenie.yaml
```         
Congratulations! Your first sample should be running through the pipeline now.  It takes right around 25 minutes for this process to complete using a single core and maxes at about 3.5 GB of memory.

We will also use `looper` to run the project pipeline locally. At the project level we can aggregate all the samples in our project (just 1 in this simple case) and view everything together.
```console
looper runp examples/meta/tutorial_refgenie.yaml
```

After the pipeline is finished, we can look through the output directory together.  We've provided a breakdown of that directory in the [browse output page](browse_output.md).

## 6. Generate an `HTML` report using `looper`

Let's take full advantage of `looper` and generate a pipeline `HTML` report that makes all our results easy to view and browse.  If you'd like to skip right to the results and see what it looks like, [check out the tutorial results](files/examples/tutorial/PEPPRO_tutorial_summary.html).  Otherwise, let's generate a report ourselves.
Using our same configuration file we used to run the samples through the pipeline, we'll now employ the `report` function of `looper`.
```console
looper report examples/meta/tutorial_refgenie.yaml
```         
That's it! Easy, right? `Looper` conveniently provides you with the location where the HTML report is produced.  You may either open the report with your preferred internet browser using the PATH provided, or we can change directories to the report's location and open it there.  Let's go ahead and change into the directory that contains the report.
```console
cd ${TUTORIAL}/processed/peppro/tutorial/
firefox PEPPRO_tutorial_summary.html
```          
The `HTML` report contains a summary page that integrates the project level summary table and any project level objects.  The status page lists all the samples in this project along with their current status, a link to their log files, the time it took to run the sample and the peak memory used during the run.  The objects page provides links to separate pages for each object type.  On each object page, all the individual samples' objects are provided.  Similarly, the samples page contains links to individual pages for each sample.  The sample pages list the individual summary statistics for that sample as well as links to log files, command logs, and summary files.  The sample pages also provide links and thumbnails for any individual objects generated for that sample.  Of course, all of these files are present in the sample directory, but the report provides easy access to them all.

# Tutorial using manually downloaded and curated genome assets

## 1: Set up folders

From an open terminal, let's first create a directory we'll use to run through this guide:
```console
mkdir peppro_tutorial
```

Let's point an environment variable to our tutorial location (change to match your local path) to link our tutorial samples with your local environment.
```console
export TUTORIAL=/path/to/your/peppro_tutorial
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
```console
git clone https://github.com/databio/peppro.git
```
Success! If you had any issues, feel free to [reach out to us with questions](contact.md).  Otherwise, let's move on to installing additional software.


## 2: Download tutorial read files

We're going to work with some files a little larger than the test data included in the pipeline so we can see all the features included in a full run of the pipeline.  Go ahead and download the [tutorial_r1.fastq.gz](http://big.databio.org/peppro/fastq/tutorial_r1.fq.gz) and [tutorial_r2.fq.gz](http://big.databio.org/peppro/fastq/tutorial_r2.fq.gz) files. 
```console
wget http://big.databio.org/peppro/fastq/tutorial_r1.fq.gz
wget http://big.databio.org/peppro/fastq/tutorial_r2.fq.gz
```

To simplify the rest of this tutorial, let's put those files in a standard location we'll use for the rest of this guide. 
```console
mv tutorial_r1.fq.gz peppro/examples/data/
mv tutorial_r2.fq.gz peppro/examples/data/
```

### 2: Get genome assets

We [recommend `refgenie` to manage all required and optional genome assets](run-bulker.md#2a-initialize-refgenie-and-download-assets). However, [`PEPPRO` can also accept file paths to any of the assets](run-bulker.md#2b-download-assets).

If you prefer not to use `refgenie`, you can also download and construct assets manually.  Because these are user-defined assets, they must exist outside of any container system. The minimum required assets for a genome includes:  
- a chromosome sizes file: a text file containing "chr" and "size" columns.  
- a [`bowtie2` genome index](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer).

Optional assets include: 

- prealignment genome assets including: `chr_sizes` and `bowtie2_index` asset for each respective genome (e.g. for `human_rDNA`)
- a TSS annotation file: a BED file containing "chr", "start", "end", "gene name", "score", and "strand" columns.
- an annotation file for pre mRNA: a BED file containing "chr", "start", "end", "gene name", "score", and "strand" columns.
- an exon annotation file: a BED file containing "chr", "start", "end", "gene name", "score", and "strand" columns.
- an intron annotation file: a BED file containing "chr", "start", "end", "gene name", "score", and "strand" columns.
- a pause index TSS annotation file: a BED file containing "chr", "start", "end", "gene name", "score", and "strand" columns.
- a pause index gene body annotation file: a BED file containing "chr", "start", "end", "gene name", "score", and "strand" columns.
- a [genomic feature annotation file](annotation.md)

You can obtain pre-built `--chrom-sizes` and `--genome-index` files from the `refgenie` servers. `Refgenie` uses algorithmically derived genome digests under-the-hood to unambiguously define genomes. That's what you'll see being used in the example below when we manually download these assets. Therefore, `2230c535660fb4774114bfa966a62f823fdb6d21acf138d4` is the digest for the human readable alias, "hg38", and `b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8` is the digest for "human_rDNA."
```console
cd ${TUTORIAL}/tools/peppro/

wget -O hg38.fasta.tgz http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/fasta?tag=default

wget -O hg38.bowtie2_index.tgz http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/bowtie2_index?tag=default

wget -O hg38.refgene_anno.tgz http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/refgene_anno?tag=default

wget -O hg38.ensembl.tgz http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/ensembl_gtf?tag=default

wget -O human_rDNA.fasta.tgz http://refgenomes.databio.org/v3/assets/archive/b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8/fasta?tag=default

wget  -O human_rDNA.bowtie2_index.tgz http://refgenomes.databio.org/v3/assets/archive/b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8/bowtie2_index?tag=default
```

Then, extract these files:
```console
tar xvf hg38.fasta.tgz
tar xvf hg38.bowtie2_index.tgz 
tar xvf hg38.refgene_anno.tgz 
tar xvf hg38.ensembl.tgz
tar xvf human_rDNA.fasta.tgz 
tar xvf human_rDNA.bowtie2_index.tgz

cd default/
wget http://big.databio.org/peppro/hg38_annotations.bed.gz
cd ../
```

## 4: Configure project files

We're going to use `looper` to analyze our data.  For that, we need to pass looper a configuration file.  This project config file describes your project. See [`looper` docs](https://looper.readthedocs.io/en/latest/) for details. A [configuration file has been provided for you in the pipeline repository itself, conveniently named `tutorial.yaml`](https://github.com/databio/peppro/blob/master/examples/tutorial/tutorial.yaml).  This configuration file also points to our sample.  In this case, we've provided a sample for you with the pipeline.  You don't have to do anything else at this point and may [skip right to running the sample if you'd like](tutorial.md#3-using-looper-to-run-the-pipeline).  Otherwise, we'll briefly touch on what those configuration files look like.
You can open the configuration file in your favorite text editor if you'd like to look closer.  For the purposes of the tutorial you may safely move past this step should you choose.
```console
nano tutorial.yaml
```
The following is what you should see in that configuration file.
```console
# Run tutorial samples through PEPPRO
name: PEPPRO_tutorial

pep_version: 2.0.0
sample_table: tutorial.csv

looper:
  output_dir: ${TUTORIAL}/processed/peppro/tutorial
  pipeline_interfaces: ${TUTORIAL}/tools/peppro/project_pipeline_interface.yaml

sample_modifiers:
  append:
    pipeline_interfaces: ${TUTORIAL}/tools/peppro/sample_pipeline_interface.yaml
  derive:
    attributes: [read1, read2]
    sources:
        R1: "${TUTORIAL}/tools/peppro/examples/data/{sample_name}_r1.fq.gz"
        R2: "${TUTORIAL}/tools/peppro/examples/data/{sample_name}_r2.fq.gz"
  imply:
    - if:
        organism: ["human", "Homo sapiens", "Human", "Homo_sapiens"]
      then:
        genome: hg38
        fasta: default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.fa
        genome_index: default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4
        chrom_sizes: default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.chrom.sizes
        prealignment_index: [human_rDNA=default/b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8]
        prealignment_names: "human_rDNA"
        TSS_name: default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_TSS.bed
        anno_name: default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.annotation.bed.gz
        pre_name: default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_pre-mRNA.bed
        exon_name: default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_exons.bed
        intron_name: default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_introns.bed
        pi_tss: default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_TSS.bed
        pi_body: default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_gene_body.bed 

```
There is also a sample annotation file referenced in our configuration file.  The sample annotation file contains metadata and other information about our sample. Just like before, [this file, named `tutorial.csv` has been provided](https://github.com/databio/peppro/blob/master/examples/meta/tutorial.csv).  You may check it out if you wish, otherwise we're all set.
If you open `tutorial.csv`, you should see the following:
```console
sample_name,protocol,organism,read1,read2,read_type
tutorial1,ATAC,human,R1,R2,paired
tutorial2,ATAC,human,R1,R2,paired
```
That's it! Let's analyze that sample!

## 5: Using `looper` to run the sample processing pipeline
Looper requires a few variables and configuration files to work for the specific user. Let's get those set up now. `Looper` uses [`divvy`](https://divvy.databio.org/) to manage computing resource configuration so that projects and pipelines can easily travel among environments. For more detailed information, [check out the `looper` docs](https://looper.readthedocs.io/en/latest/running-on-a-cluster/). Let's set it up.
```console
cd $TUTORIAL
touch compute_config.yaml
```
Open that file in your favorite text editor.  We'll add in the following example for running locally.  You'll need to edit this file further for your own setup and you can [learn more about that in the `looper` docs](https://looper.readthedocs.io/en/latest/index.html).
```console
nano compute_config.yaml
```
Paste the following into compute_config.yaml
```console
adapters:
  CODE: looper.command
  JOBNAME: looper.job_name
  CORES: compute.cores
  LOGFILE: looper.log_file
  TIME: compute.time
  MEM: compute.mem

compute_packages:
  default:
    submission_template: templates/localhost_template.sub
    submission_command: .

```
Now, let's close and save that file and create an environment variable pointing to our configuration file.
```console
export DIVCFG="${TUTORIAL}/compute_config.yaml"
```
(Remember to add `DIVCFG` to your `.bashrc` or `.profile` to ensure it persists).  
The `Looper` environment configuration file points to submission template(s) in order to know how to run a samples locally or using cluster resources.  If you'd like to learn more, check out the [`DIVCFG` configuration file and submission templates](http://code.databio.org/divvy). We're going to simply setup a local template for the purposes of this tutorial.  You can also easily create [templates for cluster or container use as well](https://github.com/pepkit/divcfg/tree/master/templates)!
Let's change to our `templates/` directory to make our first submission template.
```console
cd ${TUTORIAL}/templates/
nano localhost_template.sub
```             
Paste the following into the localhost_template.sub:
```console
#!/bin/bash

echo 'Compute node:' `hostname`
echo 'Start time:' `date +'%Y-%m-%d %T'`

{
{CODE}
} | tee {LOGFILE} --ignore-interrupts
```

Save and close that file, and return to the pipeline repository directory.
```console
cd ${TUTORIAL}/tools/peppro/
```

Now, we'll use `looper` to run the sample pipeline locally.
```console
looper run examples/meta/tutorial.yaml
```         
Congratulations! Your first sample should be running through the pipeline now.  It takes right around 25 minutes for this process to complete using a single core and maxes at about 3.5 GB of memory.

## 6: Use `looper` to run the project level pipeline
The pipeline also includes project level analyses that work on all samples concurrently.  This allows for analyses that require output produced by individual sample analysis. We'll run the project analysis much like we run the sample analysis:
```console
looper runp examples/meta/tutorial.yaml
```
This should take about a minute on the tutorial sample and will generate a `summary/` directory containing project level output in the parent project directory. 

## 7: Generate an `HTML` report using `looper`

Let's take full advantage of `looper` and generate a pipeline `HTML` report that makes all our results easy to view and browse.  If you'd like to skip right to the results and see what it looks like, [check out the tutorial results](files/examples/tutorial/PEPPRO_tutorial_summary.html).  Otherwise, let's generate a report ourselves.

Using our same configuration file we used to run the samples through the pipeline, we'll now employ the `report` function of `looper`.
```console
looper report examples/meta/tutorial.yaml
```         
That's it! Easy, right? `Looper` conveniently provides you with the location where the HTML report is produced.  You may either open the report with your preferred internet browser using the PATH returned with `looper report`, or we can change directories to the report's location and open it there.  Let's go ahead and change into the directory that contains the report.
```console
cd $TUTORIAL/processed/peppro/tutorial
firefox PEPPRO_tutorial_summary.html
```          
The `HTML` report contains a summary page that integrates the project level summary table and any project level objects including: raw aligned reads, percent aligned reads, TSS enrichment scores, and library complexity plots.  The status page lists all the samples in this project along with their current status, a link to their log files, the time it took to run the sample and the peak memory used during the run.  The objects page provides links to separate pages for each object type.  On each object page, all the individual samples' objects are provided.  Similarly, the samples page contains links to individual pages for each sample.  The sample pages list the individual summary statistics for that sample as well as links to log files, command logs, and summary files.  The sample pages also provide links and thumbnails for any individual objects generated for that sample.  Of course, all of these files are present in the sample directory, but the report provides easy access to them all.
