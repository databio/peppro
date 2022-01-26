# Run <img src="../img/peppro_logo_black.svg" alt="PEPPRO" class="img-fluid" style="max-height:35px; margin-top:-15px; margin-bottom:-10px"> with a multiple container manager.

Whether you are using `docker` or `singularity`, we have a solution to run the pipeline using containers that reduces the installation burden.

In addition to cloning the `PEPPRO` repository, this requires the installation and configuration of a single python package, our [multi-container environment manager `bulker`](https://bulker.databio.org/en/latest/). We support using `bulker` for a few reasons: 

1. It simplifies container use by wrapping the complexities of `docker` or `singularity` calls so that you can use a containerized program without even realizing you're using a container. You can call a program at the command line the same as your would *without* using bulker.
2. Similar to a dockerfile, you can distribute sets of tools *but* as a separate set of containers, not a single, unwieldy, and monolithic container.
3. Since `bulker` commands behave like native commands, a workflow becomes automatically containerized with bulker.
4. Finally, this makes bulker environments very portable, since the only requirement for native-like command use is `docker` or `singularity`.

[`Bulker` has a guide to running `PEPPRO`](https://bulker.databio.org/en/latest/peppro/), but we'll go into more detail below.

If you would still prefer using a single container, we do provide a [PEPPRO dockerfile](https://github.com/databio/peppro/blob/master/containers/peppro.Dockerfile) and support for [running the pipeline using a single, monolithic container.](run-container.md). 

## Running `PEPPRO` using `bulker`

### 1. Clone the `PEPPRO` pipeline

```console
git clone https://github.com/databio/peppro.git
```

### 2. Get genome assets

We [recommend `refgenie` to manage all required and optional genome assets](run-bulker.md#2a-initialize-refgenie-and-download-assets). However, [`PEPPRO` can also accept file paths to any of the assets](run-bulker.md#2b-download-assets).

#### 2a. Initialize `refgenie` and download assets

`PEPPRO` can utilize [`refgenie`](http://refgenie.databio.org/) assets. Because assets are user-dependent, these files must still exist outside of a container system. Therefore, we need to [install and initialize a refgenie config file.](http://refgenie.databio.org/en/latest/install/). For example:

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

`PEPPRO` also requires a `bowtie2_index` asset for any pre-alignment genomes:

```console
refgenie pull human_rDNA/fasta human_rDNA/bowtie2_index
```

#### 2b. Download assets manually

If you prefer not to use `refgenie`, you can also download and construct assets manually.  The minimum required assets for a genome includes: 
 
- a chromosome sizes file: a text file containing "chr" and "size" columns.
- a [`bowtie2` genome index](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer).

Optional assets include:  

- a TSS annotation file: a BED6 file containing "chr", "start", "end", "gene name", "score", and "strand" columns.
- a [genomic feature annotation file](annotation.md)
- a [BED6 file containing all possible TSSs used for pause-index calculation](assets.md)
- a [BED6 file containing all possible gene bodies used for pause-index calculation](assets.md)
- a [BED6 file containing gene's exon coordinates](assets.md)
- a [BED6 file containing gene's intron coordinates](assets.md)
- a [BED6 file containing premature mRNA gene coordinates](assets.md)

Even if you are *not* using `refgenie`, you can still grab premade assets for all required and optional assets from the `refgenie` servers. `Refgenie` uses algorithmically derived genome digests under-the-hood to unambiguously define genomes. That's what you'll see being used in the example below when we manually download these assets. Therefore, `2230c535660fb4774114bfa966a62f823fdb6d21acf138d4` is the digest for the human readable alias, "hg38", and `b769bcf2deaf9d061d94f2007a0e956249905c64653cb5c8` is the digest for "human_rDNA."

From the `peppro/` repository:
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
tar xvf hg38.fasta.tgz
tar xvf hg38.bowtie2_index.tgz
mv hg38.feat_annotation.gz default/2230c535660fb4774114bfa966a62f823fdb6d21acf138d4.annotation.bed.gz
tar xvf hg38.refgene_anno.tgz
tar xvf hg38.ensembl_rb.tgz
tar xvf hg38.ensembl_gtf.tgz
tar xvf human_rDNA.fasta.tgz
tar xvf human_rDNA.bowtie2_index.tgz
```

### 3. Install and configure `bulker`

Check out [the `bulker` setup guide to install bulker](https://bulker.databio.org/en/latest/install/) on your system. It is a straightforward python package with a few configuration steps required prior to use with `PEPPRO`.

### 4. Confirm installation 

After setting up your environment to run `PEPPRO` with `bulker`, you can confirm the pipeline is now executable with `bulker` using the included `checkinstall` script.  This can either be run directly from the `peppro/` repository...

```console
./checkinstall
```

or from the web:
```console
curl -sSL https://raw.githubusercontent.com/databio/peppro/checkinstall | bash
```

### 5. Load the `PEPPRO` crate

We've already produced a `bulker` crate for `PEPPRO` that requires all software needed to run the pipeline.  We can load this crate directly from the [`bulker registry`](http://hub.bulker.io/):
```console
bulker load databio/peppro:1.0.1 -r
```

### 6. Activate the `PEPPRO` crate

Now that we've loaded the `PEPPRO` crate, we need to activate that specific crate so its included tools are available.
```console
bulker activate databio/peppro:1.0.1
```
Now, you can run any of the commands in the crate as if they were natively installed, **but they're actually running in containers**!

### 7. Run the sample processing pipeline

Now we simply run the pipeline like you would with a native installation, but we wouldn't have needed to install any additional tools!

#### 7a. Run the pipeline using `looper`

Since `bulker` automatically directs any calls to required software to instead be executed in containers, we can just run our project the exact same way we would when we installed everything natively!

**Run the pipeline with looper and refgenie**
```console
looper run examples/meta/peppro_test_refgenie.yaml
```

**Run the pipeline with looper and manual asset specifications**
```console
looper run examples/meta/peppro_test.yaml
```

#### 7b. Run the pipeline at the command line

If you are using `refgenie`, but running directly at the command-line you need to specify paths to any assets that you pulled above. When [the pipeline is run with `looper`](run-bulker.md#7a-run-the-pipeline-using-looper), you can simply specify human-readable aliases to auto-populate these variables. [See the looper refgenie configuration file for an example](examples/meta/peppro_test_refgenie.yaml).

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

With a single processor, this will take around 30 minutes to complete. 

### 8. Run the project level pipeline

`PEPPRO` also includes a project-level processing pipeline to summarize the pipeline with reports on sample library complexities and count matrices across the samples in the project.

This should take < a minute on the test sample and will generate a `summary/` directory containing project level output in the parent project directory. To see more, you can [run through the extended tutorial](tutorial.md) to see this in action.

**Run the project pipeline with looper and refgenie**
```console
looper runp examples/meta/peppro_test_refgenie.yaml
```

**Run the project pipeline with looper and manual asset specifications**
```console
looper runp examples/meta/peppro_test.yaml
```

### 9. Generate an HTML report using `looper`

`Looper` can generate a pipeline HTML report that makes all our results easy to view and browse. Using the same configuration file we used to run the samples through the pipeline, we'll now employ the `report` function of `looper`.

**Generate the HTML report with looper and refgenie managed assets**
```console
looper report examples/meta/peppro_test_refgenie.yaml
```

**Generate the HTML report with looper and manual asset specifications**
```console
looper report examples/meta/peppro_test.yaml
```

