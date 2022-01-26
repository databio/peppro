# Genome assets

`PEPPRO` can use either manually constructed or `refgenie` managed assets. `Refgenie` streamlines sample processing, where once assets are built by `refgenie` there is minimal argument calls to `PEPPRO` to use all assets. Pipeline assets include:  

**Required**  

| `PEPPRO` argument | `refgenie` asset name                                                                                                                                      | Description                                                                                           |
|--------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------|
| `--genome-index`   | [`bowtie2_index`](http://refgenie.databio.org/en/latest/available_assets/#bowtie2_index)                                                                   | A genome index file constructed from `bowtie2-build`                                                  |
| `--chrom-sizes`    | With `refgenie`, this asset is built automatically when you build/pull the [`fasta`](http://refgenie.databio.org/en/latest/available_assets/#fasta) asset. | A text file containing "chr" and "size" columns.                                                      |

**Optional**  

| `PEPPRO` argument     | `refgenie` asset name                                                                                                                                    | Description                                                                                                                                                                             |
|------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `--prealignment-names` | Human readable genome alias(es) for `refgenie` managed `bowtie2_index` asset(s).                                                                         | A space-delimited list of genome names. *e.g.* ["rCRSd", "human_repeats"]                                                                                                               |
| `--prealignment-index` | [`bowtie2_index`](http://refgenie.databio.org/en/latest/available_assets/#bowtie2_index)                                                                 | A genome index file constructed from `bowtie2-build`. Used for manually pointing to prealignment genome indices when using `bowtie2` (default) for alignment.                           |
| `--TSS-name`           | [`refgene_anno`](http://refgenie.databio.org/en/latest/available_assets/#refgene_anno). `refgenie` `build/pull` the TSS annotation file with this asset. | Transcription start site (TSS) annotations. *e.g.* [refGene.txt.gz](https://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz)                                            |
| `--anno-name`          | [`feat_annotation`](annotation.md)                                                                                                                       | A BED-style file with "chr", "start", "end", "genomic feature name", "score" and "strand" columns.                                                                                      |
| `--pi-tss`             | [`ensembl_gtf.ensembl_tss`](http://refgenie.databio.org/en/latest/available_assets/#ensembl_gtf)                                                         | A derived asset from an Ensembl GTF file. Represents all possible TSSs.                                                                                                                 |
| `--pi-body`            | [`ensembl_gtf.ensembl_gene_body`](http://refgenie.databio.org/en/latest/available_assets/#ensembl_gtf)                                                   | A derived asset from an Ensembl GTF file. Represents all possible gene body coordinates.                                                                                                |
| `--pre-name`           | [`refgene_anno.refgene_pre_mRNA`](http://refgenie.databio.org/en/latest/available_assets/#refgene_anno)                                                  | Asset derived from a refGene annotation file. Represents premature mRNA coordinates.                                                                                                    |
| `--exon-name`          | [`refgene_anno.refgene_exon`](http://refgenie.databio.org/en/latest/available_assets/#refgene_anno)                                                      | Asset derived from a refGene annotation file. Represents all exon coordinates.                                                                                                          |
| `--intron-name`        | [`refgene_anno.refgene_intron`](http://refgenie.databio.org/en/latest/available_assets/#refgene_anno)                                                    | Asset derived from a refGene annotation file. Represents all intron coordinates.                                                                                                        |
| `--fasta`              | [`fasta`](https://refgenie.databio.org/en/latest/available_assets/#fasta) The `fasta` asset.                                                             | A genome fasta file. Required for `--sob` argument.                                                                                                                                     |
| `--search-file`        | [`tallymer_index`](https://refgenie.databio.org/en/latest/available_assets/#tallymer_index) The `search_file` is built from this `refgenie` asset.       | File used to search an index of k-mers in the genome of the same size as input read lengths. Only required for `--sob` argument                                                         |

## Using `refgenie` managed assets

`PEPPRO` can utilize [`refgenie`](http://refgenie.databio.org/) assets. Because assets are user-dependent, these files must be available natively. Therefore, you need to [install and initialize a refgenie config file.](http://refgenie.databio.org/en/latest/install/). For example:

```console
pip install refgenie
export REFGENIE=/path/to/your_genome_folder/genome_config.yaml
refgenie init -c $REFGENIE
```

Add the `export REFGENIE` line to your `.bashrc` or `.profile` to ensure it persists. 

Next, pull the assets you need. Replace `hg38` in the example below if you need to use a different genome assembly. If these assets are not available automatically for your genome of interest, then you'll need to [build them](annotation.md). Download all standard assets for `hg38` like so:

```console
refgenie pull hg38/fasta hg38/bowtie2_index hg38/refgene_anno hg38/ensembl_gtf hg38/ensembl_rb
refgenie build hg38/feat_annotation
```

`PEPPRO` also requires a `fasta` and `bowtie2_index` asset for any prealignment genomes:

```console
refgenie pull human_rDNA/fasta human_rDNA/bowtie2_index
```

Furthermore, you can [learn more about using `seqOutBias` and the required `tallymer_index` here](sob.md).

### Example using `refgenie` managed assets

When using `refgenie`, you only need to provide the `--genome` and `--prealignment-names` argument to provide the pipeline with every required index and optional annotation file that exists for those genomes. This means, the TSS file, feature annotation file, and blacklist will all be used without needing to directly specify the paths to these files.

From the `peppro/` repository directory:
```console
looper run examples/meta/peppro_test_refgenie.yaml
```

## Using manually managed assets

Assets may also be managed manually and specified directly to the pipeline.  While this frees you from needing `refgenie` installed and initialized, it does require a few more arguments to be specified.

The TSS annotation file may be specified using `--TSS-name </path/to/your_TSS_annotations.bed>`. This file is a `BED6` (e.g. chr, start, end, name, score, strand) formatted file.

The `feat_annotation` asset may also be directly specified using `--anno-name </path/to/your_custom_feature_annotations.bed.gz>`.  Read [more about using custom reference data](annotation.md).

The `pi_tss` asset, representing all possible TSSs for calculating the pause index, may be directly specified using `--pi-tss`. This file is a `BED6` (e.g. chr, start, end, name, score, strand) formatted file.

The `pi_body` asset, representing all possible gene bodies for calculating the pause index, may be directly specified using `--pi-body`. This file is a `BED6` (e.g. chr, start, end, name, score, strand) formatted file.

The `pre_name` asset, representing premature mRNA sequence coordinates, may be directly specified using `--pre-name`. This file is a `BED6` (e.g. chr, start, end, name (a gene name), score, strand) formatted file.

The `exon_name` asset, representing gene exon coordinates, may be directly specified using `--exon-name`. This file is a `BED6` (e.g. chr, start, end, name (the name of the gene the exon is from), score, strand) formatted file.

The `intron_name` asset, representing gene intron coordinates, may be directly specified using `--intron-name`. This file is a `BED6` (e.g. chr, start, end, name (the name of the gene the intron is from), score, strand) formatted file.

### Example using manually managed assets

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

Then, extract these files to the `peppro/` parent directory:
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

From the `peppro/` repository folder (using the manually downloaded genome assets):
```console
looper run examples/meta/peppro_test.yaml
```
