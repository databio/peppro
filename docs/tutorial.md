# <img src="../img/peppro_logo.svg" alt="PEPPRO" class="img-fluid" style="max-height:50px; margin-top:-15px; margin-bottom:-10px"> pipeline step-by-step guide

In this guide, we'll walk you through the step by step procedure of running a tutorial PRO-seq dataset through the pipeline.  The output from this process is the same as you see in the [example PRO-seq output](browse_output.md) we've provided.  To use this tutorial, you should have a basic familiarity with [working in a command line driven environment](http://matt.might.net/articles/basic-unix/). You also need to have already installed `PEPPRO` prerequisites, which you can do following the [basic installation instructions](install.md).

## 1: Set up folders 

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

## 2: Download tutorial read files

We're going to work with some files a little larger than the test data included in the pipeline so we can see all the features included in a full run of the pipeline.  Go ahead and download the [tutorial_r1.fastq.gz](http://big.databio.org/peppro/tutorial_r1.fq.gz) and [tutorial_r2.fq.gz](http://big.databio.org/peppro/tutorial_r2.fastq.gz) files. 
```console
wget http://big.databio.org/peppro/tutorial_r1.fq.gz
wget http://big.databio.org/peppro/tutorial_r2.fq.gz
```

To simplify the rest of this tutorial, let's put those files in a standard location we'll use for the rest of this guide. 
```console
mv tutorial_r1.fq.gz peppro/examples/data/
mv tutorial_r2.fq.gz peppro/examples/data/
```

## 3: Configure project files

We're going to use `looper` to analyze our data.  For that, we need to pass looper a configuration file.  This project config file describes your project. See [`looper` docs](https://looper.readthedocs.io/en/latest/) for details. A configuration file has been provided for you in the pipeline itself, conveniently named `tutorial.yaml`.  This configuration file also points to our sample.  In this case, we've provided a sample for you with the pipeline.  You don't have to do anything else at this point and may [skip right to running the sample if you'd like](tutorial.md#3-using-looper-to-run-the-pipeline).  Otherwise, we'll briefly touch on what those configuration files look like.

You can open the configuration file in your favorite text editor if you'd like to look closer.  For the purposes of the tutorial you may safely move past this step should you choose.
```
nano tutorial.yaml
```
The following is what you should see in that configuration file.
```
# Run tutorial sample through PEPPRO
name: tutorial

metadata:
  sample_annotation: "peppro_tutorial.csv"
  output_dir: "$PROCESSED/peppro_tutorial"
  pipeline_interfaces: "$CODE/peppro/pipeline_interface.yaml" 
  
derived_columns: [read1, read2]

data_sources:
  R1: "$CODE/peppro/examples/data/{sample_name}_r1.fq.gz"
  R2: "$CODE/peppro/examples/data/{sample_name}_r2.fq.gz"
  
implied_columns:
  organism:
    human:
      genome: hg38
      prealignments: human_rDNA rCRSd
```
There is also a sample annotation file referenced in our configuration file.  The sample annotation file contains metadata and other information about our sample. Just like before, this file, named `tutorial.csv` has been provided.  You may check it out if you wish, otherwise we're all set.
If you open `tutorial.csv`, you should see the following:
```
sample_name,organism,protocol,read_type,read1,read2
tutorial,human,PROSEQ,paired,R1,R2
```
That's it! Let's analyze that sample!


## 4: Create environment variables

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

## 5: Use `looper` to run the pipeline
Looper requires a few variables and configuration files to work for the specific user. Let's get those set up now. `Looper` uses [`divvy`](https://divvy.databio.org/) to manage computing resource configuration so that projects and pipelines can easily travel among environments. For more detailed information, [check out the `looper` docs](https://looper.readthedocs.io/en/latest/cluster-computing/). Let's set it up.
```
cd /path/to/peppro_tutorial/
touch compute_config.yaml
```
Open that file in your favorite text editor.  We'll add in the following example for running locally.  You'll need to edit this file further for your own setup and you can [learn more about that in the `looper` docs](https://looper.readthedocs.io/en/latest/index.html).
```
nano compute_config.yaml
```
Paste the following into compute_config.yaml
```
compute:
  default:
    submission_template: templates/localhost_template.sub
    submission_command: sh
```
Now, let's close and save that file and create an environment variable pointing to our configuration file.
```
export DIVCFG="/path/to/peppro_tutorial/compute_config.yaml"
```
(Remember to add `DIVCFG` to your `.bashrc` or `.profile` to ensure it persists).
The `looper` environment configuration file points to submission template(s) in order to know how to run a samples locally or using cluster resources.  If you'd like to learn more, check out the [`DIVCFG` configuration file and submission templates](https://divvy.databio.org/). We're going to simply setup a local template for the purposes of this tutorial.  You can also easily create [templates for cluster or container use as well](https://github.com/pepkit/divcfg/tree/master/templates)!
Let's change to our `templates/` directory to make our first submission template.
```
cd /path/to/peppro_tutorial/templates/
nano localhost_template.sub
```             
Paste the following into the localhost_template.sub:
```
#!/bin/bash

echo 'Compute node:' `hostname`
echo 'Start time:' `date +'%Y-%m-%d %T'`

{
{CODE}
} | tee {LOGFILE} --ignore-interrupts
```

Save and close that file, and return to our main tutorial directory.
```
cd ../
```
Now, we'll use `looper` to run the sample locally.
```
looper run tutorial.yaml
```         
Congratulations! Your first sample should be running through the pipeline now.  It takes right around 25 minutes for this process to complete using a single core and maxes at about 3.5 GB of memory.

After the pipeline is finished, we can look through the output directory together.  We've provided a breakdown of that directory in the [browse output page](/browse_output/).


## 6: Generate an `HTML` report using `looper`

Let's take full advantage of `looper` and generate a pipeline `HTML` report that makes all our results easy to view and browse.  If you'd like to skip right to the results and see what it looks like, [check out the tutorial results](../files/examples/tutorial/tutorial_summary.html).  Otherwise, let's generate a report ourselves.
Using our same configuration file we used to run the samples through the pipeline, we'll now employ the `summarize` function of `looper`.
```
looper summarize tutorial.yaml
```         
That's it! Easy, right? `Looper` conveniently provides you with the location where the HTML report is produced.  You may either open the report with your preferred internet browser using the PATH provided, or we can change directories to the report's location and open it there.  Let's go ahead and change into the directory that contains the report.
```
cd /path/to/peppro_tutorial/processed/tutorial/
firefox tutorial_summary.html
```          
The `HTML` report contains a summary page that integrates the project level summary table and any project level objects including: raw aligned reads, percent aligned reads, and TSS enrichment scores.  The status page lists all the samples in this project along with their current status, a link to their log files, the time it took to run the sample and the peak memory used during the run.  The objects page provides links to separate pages for each object type.  On each object page, all the individual samples' objects are provided.  Similarly, the samples page contains links to individual pages for each sample.  The sample pages list the individual summary statistics for that sample as well as links to log files, command logs, and summary files.  The sample pages also provide links and thumbnails for any individual objects generated for that sample.  Of course, all of these files are present in the sample directory, but the report provides easy access to them all.