# Run <img src="../img/peppro_logo_black.svg" alt="PEPPRO" class="img-fluid" style="max-height:35px; margin-top:-15px; margin-bottom:-10px"> in a container.

A popular approach is installing all dependencies in a container and just use that single container. This container can be used with either `docker` or `singularity`. You can run `PEPPRO` as an individual pipeline on a single sample using the container with `docker run` or `singularity exec`. Or, you can rely on `looper`, which is already set up to run any pipeline in existing containers using the `divvy` templating system. 

## Running `PEPPRO` using a single, monolithic container.

### 1: Clone the `PEPPRO` pipeline

```console
git clone https://github.com/databio/peppro.git
```

### 2: Get genome assets

We [recommend `refgenie` to manage all required and optional genome assets](run-container.md#2a-initialize-refgenie-and-download-assets). However, [`PEPPRO` can also accept file paths to any of the assets](run-container.md#2b-download-assets).

#### 2a: Initialize `refgenie` and download assets

`PEPPRO` can use [`refgenie`](http://refgenie.databio.org/) assets for alignment and annotation. Because assets are user-dependent, these files must still exist outside of a container system. We need to [install and initialize a refgenie config file.](http://refgenie.databio.org/en/latest/install/). For example:

```console
pip install refgenie
export REFGENIE=/path/to/your_genome_folder/genome_config.yaml
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

#### 2b: Download assets manually

If you prefer not to use `refgenie`, you can also download and construct assets manually.  Again, because these are user-defined assets, they must exist outside of any container system. The minimum required assets for a genome includes:  
- a chromosome sizes file: a text file containing "chr" and "size" columns.  
- a [`bowtie2` genome index](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer).
- an [ensembl_gtf](http://refgenie.databio.org/en/latest/available_assets/#ensembl_gtf) asset used to build other derived assets including a comprehensive TSS annotation and gene body annotation.
- an [ensembl_rb] (http://refgenie.databio.org/en/latest/available_assets/#ensembl_rb) asset containing known genomic features such as promoters and used to produce derived assets such as genomic feature annotations.
- a [refgene_anno](http://refgenie.databio.org/en/latest/available_assets/#refgene_anno) asset used to produce derived assets including transcription start sites (TSSs), exons, introns, and premature mRNA sequences.
- a [genomic feature annotation file](annotation.md) (which may also be built locally through the `refgenie build <genome_name>/feat_annotation`)

You can still obtain the pre-constructed assets from the [`refgenie` servers](http://refgenomes.databio.org/v3/genomes/splash/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4). `Refgenie` uses algorithmically derived genome digests under-the-hood to unambiguously define genomes. That's what you'll see being used in the example below when we manually download these assets. Therefore, `2230c535660fb4774114bfa966a62f823fdb6d21acf138d4` is the digest for the human readable alias, "hg38", and `b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8` is the digest for "human_rDNA."
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

Then, extract these files:
```console
tar xf hg38.fasta.tgz
tar xf hg38.bowtie2_index.tgz
tar xf hg38.ensembl_gtf.tgz
tar xf hg38.ensembl_rb.tgz
tar xf hg38.refgene_anno.tgz
mv hg38.feat_annotation.gz default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.annotation.bed.gz
tar xf human_rDNA.fasta.tgz
tar xf human_rDNA.bowtie2_index.tgz
```

### 3. Pull the container image.

**Docker**: You can pull the docker [databio/peppro image](https://hub.docker.com/r/databio/peppro/) from `dockerhub` like this:

```console
docker pull databio/peppro
```

Or build the image using the included `Dockerfile` (you can use a recipe in the included `Makefile` in the `peppro/` repository):
```console
make docker
```

**Singularity**: You can [download the `singularity` image](http://big.databio.org/simages/peppro) or build it from the docker image using the `Makefile`:
```console
make singularity
```

Now you'll need to tell the pipeline where you saved the singularity image. You can either create an environment variable called `$SIMAGES` that points to the folder where your image is stored, or you can tweak the `pipeline_interface.yaml` file so that the `compute.singularity_image` attribute is pointing to the right location on disk.

### 6. Confirm installation 

After setting up your environment to run `PEPPRO` using containers, you can confirm the pipeline is now executable with your container system using the included `checkinstall` script.  This can either be run directly from the `peppro/` repository...

```console
./checkinstall
```

or from the web:
```console
curl -sSL https://raw.githubusercontent.com/databio/peppro/checkinstall | bash
```

### 4. Run individual samples in a container

Individual jobs can be run in a container by simply running the `peppro.py` command through `docker run` or `singularity exec`. You can run containers either on your local computer, or in an HPC environment, as long as you have `docker` or `singularity` installed. You will need to include any volumes that contain data required by the pipeline. For example, to utilize `refgenie` assets you'll need to ensure the volume containing those files is available. In the following example, we are including an environment variable (`$GENOMES`) which points to such a directory.

For example, run it locally in `singularity` like this:
```console
singularity exec $SIMAGES/peppro pipelines/peppro.py --help
```

With `docker`, you can use:
```console
docker run --rm -it databio/peppro pipelines/peppro.py --help
```

### 5. Running multiple samples in a container with looper

To run multiple samples in a container, you simply need to configure `looper` to use a container-compatible template. The looper documentation has instructions for [running jobs in containers](http://looper.databio.org/en/latest/containers/).

### Container details 

#### Using `docker`
The pipeline has been successfully run in both a `Linux` and `MacOS` environment. With `docker` you need to bind mount your volume that contains the pipeline and your genome assets locations, as well as provide the container the same environment variables your host environment is using.

In the first example, we're mounting our home user directory (`/home/jps3ag/`) which contains the parent directories to our genome assets and to the pipeline itself. We'll also provide the pipeline environment variables, such as `$HOME`.

Here's that example command in a Linux environment to run the test example through the pipeline (using the manually downloaded genome assets):
```console
docker run --rm -it --volume /home/jps3ag/:/home/jps3ag/ \
  -e HOME='/home/jps3ag/' \
  databio/peppro \
  /home/jps3ag/src/peppro/pipelines/peppro.py --single-or-paired single \
  --prealignment-index human_rDNA=default/b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8 \
  --genome hg38 \
  --genome-index /home/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4 \
  --chrom-sizes /home/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.chrom.sizes \
  --sample-name test \
  --input /home/jps3ag/src/peppro/examples/data/test_r1.fq.gz \
  --TSS-name /home/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_TSS.bed \
  --anno-name /home/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.annotation.bed.gz \
  --pre-name /home/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_pre-mRNA.bed \
  --exon-name /home/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_exons.bed \
  --intron-name /home/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_introns.bed \
  --pi-tss /home/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_TSS.bed \
  --pi-body /home/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_gene_body.bed \
  -O $HOME/peppro_test
```

In this second example, we'll perform the same command in a `MacOS` environment using [`Docker` for `Mac`](https://docs.docker.com/desktop/mac/install/). 

This necessitates a few minor changes to run that same example:

- replace `/home/` with `/Users/` format
- e.g. `--volume /Users/jps3ag/:/Users/jps3ag/`

Be sure to [allocate sufficient memory](https://docs.docker.com/desktop/mac/#resources) (6-8GB should generally be adequate) in Docker for Mac.

```console
docker run --rm -it --volume /Users/jps3ag/:/Users/jps3ag/ \
  -e HOME="/Users/jps3ag/" \
  databio/peppro \
  /Users/jps3ag/src/peppro/pipelines/peppro.py --single-or-paired single \
  --prealignment-index human_rDNA=/Users/jps3ag/src/peppro/default/b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8 \
  --genome hg38 \
  --genome-index /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4 \
  --chrom-sizes /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.chrom.sizes \
  --sample-name test \
  --input /Users/jps3ag/src/peppro/examples/data/test_r1.fq.gz \
  --TSS-name /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_TSS.bed \
  --anno-name /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.annotation.bed.gz \
  --pre-name /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_pre-mRNA.bed \
  --exon-name /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_exons.bed \
  --intron-name /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_introns.bed \
  --pi-tss /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_TSS.bed \
  --pi-body /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_gene_body.bed \
  -O peppro_test
```

#### Using `singularity`

First, build a singularity container from the docker image and create a running instance:
```console
singularity build peppro docker://databio/peppro:latest
singularity instance start -B /home/jps3ag/:/home/jps3aq/ peppro peppro_instance
```

Second, run your command.
```console
singularity exec instance://peppro_instance \
  /home/jps3ag/src/peppro/pipelines/peppro.py --single-or-paired single \
  --prealignment-index human_rDNA=/Users/jps3ag/src/peppro/default/b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8 \
  --genome hg38 \
  --genome-index /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4 \
  --chrom-sizes /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.chrom.sizes \
  --sample-name test \
  --input /home/jps3ag/src/peppro/examples/data/test_r1.fq.gz \
  --TSS-name /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_TSS.bed \
  --anno-name /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.annotation.bed.gz \
  --pre-name /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_pre-mRNA.bed \
  --exon-name /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_exons.bed \
  --intron-name /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_introns.bed \
  --pi-tss /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_TSS.bed \
  --pi-body /Users/jps3ag/src/peppro/default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4_ensembl_gene_body.bed \
  -O peppro_test
```

Third, close your instance when finished.
```
singularity instance stop peppro_instance
```

