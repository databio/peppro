# Run <img src="../img/peppro_logo_black.svg" alt="PEPPRO" class="img-fluid" style="max-height:35px; margin-top:-15px; margin-bottom:-10px"> in a conda environment.

We also enable setup of the pipeline using `conda`. As with container-based approaches, some native installation is required for complete setup. 

## 1. Clone the `PEPPRO` pipeline

```console
git clone https://github.com/databio/peppro.git
```

## 2. Install bioinformatic tools

Be prepared for this initial installation process to take more than an hour to complete.

From the `peppro/` repository directory:
```{bash}
conda env create -f requirements-conda.yml
```

Note: The subsequent steps all assume you have installed using `conda`.  Alternatively, you can [follow instructions to install each individual program natively](detailed-install.md).

## 3. Install python packages

`PEPPRO` uses several Python packages under the hood. Not all of these are available through `conda`, so we'll ensure they are installed ourselves to the `peppro` `conda` environment. From the `peppro/` directory:

```{bash}
conda activate peppro
unset PYTHONPATH
python -m pip install --ignore-installed --upgrade -r requirements.txt
```

## 4. Install R packages

`PEPPRO` uses `R` to generate quality control and read/peak annotation plots. We have packaged the `peppro` specific `R` code into a supporting package called [PEPPROr](https://github.com/databio/peppro/tree/master/PEPPROr). The `PEPPROr` package relies on a few additional packages which can be installed to the `conda` environment.

To ensure these packages are installed to the `peppro` `conda` environment, make sure to point your `R_LIBS` environment variable to the `conda` environment `R` library. For example:
```{bash}
conda activate peppro
unset R_LIBS
export R_LIBS="$CONDA_PREFIX/lib/R/library"
```

From the `peppro/` directory, open `R` and install the following packages:
```{R}
install.packages("optigrab")
devtools::install_github("databio/GenomicDistributions")
install.packages("http://big.databio.org/GenomicDistributionsData/GenomicDistributionsData_0.0.2.tar.gz", repos=NULL)
devtools::install(file.path("PEPPROr/"), dependencies=TRUE, repos="https://cloud.r-project.org/")
```

## 5. Get genome assets

### 5a. Initialize `refgenie` and download assets

`PEPPRO` can utilize [`refgenie`](http://refgenie.databio.org/) assets. Because assets are user-dependent, these files must still be available natively. Therefore, we need to [install and initialize a refgenie config file.](http://refgenie.databio.org/en/latest/install/). For example:

```console
pip install refgenie
export REFGENIE=genome_config.yaml
refgenie init -c $REFGENIE
```

Add the `export REFGENIE` line to your `.bashrc` or `.profile` to ensure it persists. 

Next, pull the assets you need. Replace `hg38` in the example below if you need to use a different genome assembly. If these assets are not available automatically for your genome of interest, then you'll need to [build them](annotation.md). Download these required assets with this command:

```console
refgenie pull hg38/fasta hg38/bowtie2_index hg38/refgene_anno hg38/ensembl_gtf hg38/ensembl_rb
refgenie build hg38/feat_annotation
```

`PEPPRO` also requires a `fasta` and `bowtie2_index` asset for any pre-alignment genomes:

```console
refgenie pull human_rDNA/fasta human_rDNA/bowtie2_index
```

### 5b. Download assets manually

If you prefer not to use `refgenie`, you can also download [assets](assets.md) manually. To realize the full potential of the pipeline, you will need the following:
 
 - a chromosome sizes file: a text file containing "chr" and "size" columns.  
 - a [`bowtie2` genome index](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer).
 - an [ensembl_gtf](http://refgenie.databio.org/en/latest/available_assets/#ensembl_gtf) asset used to build other derived assets including a comprehensive TSS annotation and gene body annotation.
 - an [ensembl_rb] (http://refgenie.databio.org/en/latest/available_assets/#ensembl_rb) asset containing known genomic features such as promoters and used to produce derived assets such as genomic feature annotations.
 - a [refgene_anno](http://refgenie.databio.org/en/latest/available_assets/#refgene_anno) asset used to produce derived assets including transcription start sites (TSSs), exons, introns, and premature mRNA sequences.
 - a [genomic feature annotation file](annotation.md)
 
Even if you are *not* using `refgenie`, you can still grab these assets for all required and optional assets from the `refgenie` servers. `Refgenie` uses algorithmically derived genome digests under-the-hood to unambiguously define genomes. That's what you'll see being used in the example below when we manually download these assets. Therefore, `2230c535660fb4774114bfa966a62f823fdb6d21acf138d4` is the digest for the human readable alias, "hg38", and `b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8` is the digest for "human_rDNA."

From within the `peppro/` repository:
```console
wget -O hg38.fasta.tgz http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/fasta?tag=default
wget  -O hg38.bowtie2_index.tgz http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/bowtie2_index?tag=default
wget  -O hg38.ensembl_gtf.tgz http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/ensembl_gtf?tag=default
wget  -O hg38.ensembl_rb.tgz http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/ensembl_rb?tag=default
wget  -O hg38.refgene_anno.tgz http://refgenomes.databio.org/v3/assets/archive/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4/refgene_anno?tag=default
wget -O hg38.feat_annotation.gz http://big.databio.org/peppro/hg38_annotations.bed.gz
wget  -O human_rDNA.fasta.tgz http://refgenomes.databio.org/v3/assets/archive/b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8/fasta?tag=default
wget  -O human_rDNA.bowtie2_index.tgz http://refgenomes.databio.org/v3/assets/archive/b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8/bowtie2_index?tag=default
```

Then, extract those files:
```console
tar xvf hg38.fasta.tgz
tar xvf hg38.bowtie2_index.tgz
mv hg38.feat_annotation.gz default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.annotation.bed.gz
tar xvf hg38.refgene_anno.tgz
tar xvf hg38.ensembl_rb.tgz
tar xvf hg38.ensembl_gtf.tgz
tar xvf human_rDNA.fasta.tgz
tar xvf human_rDNA.bowtie2_index.tgz
```

## 6. Confirm installation 

After setting up your environment to run `PEPPRO` with `conda`, you can confirm the pipeline is executable with `conda` using the included `checkinstall` script.  This can either be run directly from the `peppro/` repository...

```console
./checkinstall
```

or from the web:
```console
curl -sSL https://raw.githubusercontent.com/databio/peppro/checkinstall | bash
```

## 7. Run the sample processing pipeline

Now we can run the pipeline in the `peppro` conda environment. The easiest approach is to use `looper`, but you can also run the pipeline for a single sample directly at the command line.

### 7a. Run the pipeline using `looper`

`PEPPRO` can utilize a [pipeline submission engine called `looper`](http://looper.databio.org/en/latest/) to run the pipeline across each sample in a project. We can use the `-d` argument to first try a dry run, which will create job scripts for every sample in a project, but will not execute them.

**Run the pipeline with looper and refgenie**
```console
looper run examples/meta/peppro_test_refgenie.yaml
```

**Run the pipeline with looper and manual asset specifications**
```console
looper run examples/meta/peppro_test.yaml
```

There are lots of other cool things you can do with `looper`, like the dry runs, or report results, check on pipeline run status, clean intermediate files to save disk space, lump multiple samples into one job, and more. For details, consult the [looper docs](http://looper.databio.org/).

### 7b. Run the pipeline at the command line

If you are using `refgenie`, but running directly at the command-line you need to specify paths to any assets that you pulled above. When [the pipeline is run with `looper`](run-conda.md#7a-run-the-pipeline-using-looper), you can simply specify human-readable aliases to auto-populate these variables. [See the looper refgenie configuration file for an example](examples/meta/peppro_test_refgenie.yaml).

You can grab the path to the minimally required `--chrom-sizes` and `--genome-index` files as follows:
```console
refgenie seek hg38/fasta.chrom_sizes
refgenie seek hg38/bowtie2_index.dir
```

And if you are using pre-alignments, you need the genome index for any pre-alignment genomes, `--prealignment-index`:
```console
refgenie seek human_rDNA/bowtie2_index.dir
```

For the full potential of the pipeline, you'll also need the file paths for the following assets:

| pipeline argument | `refgenie` command to retrieve file path           |
|-------------------|----------------------------------------------------|
| `--TSS-name`      | `refgenie seek hg38/refgene_anno.refgene_tss`      |
| `--anno-name`     | `refgenie seek hg38/feat_annotation`               |
| `--pre-name`      | `refgenie seek hg38/refgene_anno.refgene_pre_mRNA` |
| `--exon-name`     | `refgenie seek hg38/refgene_anno.refgene_exon`     |
| `--intron-name`   | `refgenie seek hg38/refgene_anno.refgene_intron`   |
| `--pi-tss`        | `refgenie seek hg38/ensembl_gtf.ensembl_tss`       |
| `--pi-body`       | `refgenie seek hg38/ensembl_gtf.ensembl_gene_body` |

You'll need to update the paths to the assets to reflect the results from `refgenie seek`. Below is an example where all those assets are local to the `peppro/` repository.

From the `peppro/` repository folder (using `refgenie` managed genome assets file paths):
```console
pipelines/peppro.py --single-or-paired single \
  --prealignment-index human_rDNA=default/b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8 \
  --genome-index default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4 \
  --chrom-sizes default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.chrom.sizes \
  --genome hg38 \
  --sample-name test \
  --input examples/data/test_r1.fq.gz \
  --TSS-name default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_TSS.bed \
  --anno-name default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.annotation.bed.gz \
  --pre-name default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_pre-mRNA.bed \
  --exon-name default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_exons.bed \
  --intron-name default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_introns.bed \
  --pi-tss default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_TSS.bed \
  --pi-body default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_gene_body.bed \
  -O peppro_test
```

In the previous example, we used `refgenie` assets that we placed in the same location as if we manually downloaded assets to the `peppro/` repository, so the file paths here look the same.  From the `peppro/` repository folder (using the manually downloaded genome assets):
```console
pipelines/peppro.py --single-or-paired single \
  --prealignment-index human_rDNA=default/b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8 \
  --genome hg38 \
  --genome-index default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4 \
  --chrom-sizes default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.chrom.sizes \
  --sample-name test \
  --input examples/data/test_r1.fq.gz \
  --TSS-name default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_TSS.bed \
  --anno-name default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.annotation.bed.gz \
  --pre-name default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_pre-mRNA.bed \
  --exon-name default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_exons.bed \
  --intron-name default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_introns.bed \
  --pi-tss default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_TSS.bed \
  --pi-body default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_gene_body.bed \
  -O peppro_test
```

## 8. Use `looper` to run the project level pipeline

`PEPPRO` also includes a project-level processing pipeline to summarize the pipeline with reports on sample library complexities and count matrices across the samples in the project.

**Run the project pipeline with looper and refgenie managed assets**
```console
looper runp examples/meta/peppro_test_refgenie.yaml
```

**Run the project pipeline with looper and manual asset specifications**
```console
looper runp examples/meta/peppro_test.yaml
```

This should take < a minute on the test sample and will generate a `summary/` directory containing project level output in the parent project directory. In this small example, there won't be a consensus peak set or count table because it is only a single sample. To see more, you can [run through the extended tutorial](tutorial.md) to see this in action.

## 9. Generate an HTML report using `looper`

`Looper` can generate a pipeline HTML report that makes all our results easy to view and browse. Using the same configuration file we used to run the samples through the pipeline, we'll now employ the `report` function of `looper`.

**Generate the HTML report with looper and refgenie managed assets**
```console
looper report examples/meta/peppro_test_refgenie.yaml
```

**Generate the HTML report with looper and manual asset specifications**
```console
looper report examples/meta/peppro_test.yaml
```