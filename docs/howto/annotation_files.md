# Download or create annotation files for <img src="../../img/peppro_logo_black.svg" alt="PEPPRO" class="img-fluid" style="max-height:35px; margin-top:-15px; margin-bottom:-10px"> 


For each annotation type (TSS, CpA sites, premature mRNA, or general features), we provide [downloadable defaults](http://big.databio.org/peppro/) for common genomes.  You may also recreate these yourself as described below.

### TSS

To calculate [TSS enrichments](../glossary.md), you will need a [TSS annotation file](http://big.databio.org/refgenomes/) in your reference genome directory.  If a pre-built version for your genome of interest isn't present, you can quickly create that file yourself. In the reference genome directory, you can perform the following commands for in this example, `hg38`:
```console
wget -O hg38_TSS_full.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz \
zcat hg38_TSS_full.txt.gz | \
  awk  '{if($4=="+"){print $3"\t"$5"\t"$5"\t"$4"\t"$13}else{print $3"\t"$6"\t"$6"\t"$4"\t"$13}}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_TSS.tsv
```
This asset (`tss_annotation`) needs to be [included in your `$REFGENIE` configuration file](#Example_PEPPRO_REFGENIE_configuration_file) for the pipeline to detect it automatically.  Alternatively, you can use the `--TSS-name` pipeline option to provide a path directly to this file.

### Pause index annotation (PI)

To calculate [pause indicies](../glossary.md), you will need a [PI annotation file](http://big.databio.org/refgenomes/) in your reference genome directory.  If a pre-built version for your genome of interest isn't present, you can quickly create that file yourself. In the reference genome directory, you can perform the following commands for in this example, `hg38`:
```console
wget -O hg38_TSS_full.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz \
zcat hg38_TSS_full.txt.gz | \
  awk  '{if($4=="+"){print $3"\t"$5-80"\t"$5-20"\t"$13"\t"1"\t"$4"\t"$3"\t"$5-500"\t"$6"\t"$13"\t"1"\t"$4}else{print $3"\t"$6+20"\t"$6+80"\t"$13"\t"1"\t"$4"\t"$3"\t"$5"\t"$6+500"\t"$13"\t"1"\t"$4}}' | \
  awk -v OFS='\t' '$2<0{$2=0}1' | \
  awk -v OFS='\t' '$3<0{$3=0}1' | \
  awk -v OFS='\t' '$7<0{$7=0}1' | \
  awk -v OFS='\t' '$8<0{$8=0}1' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_PI.tsv
```
This asset (`pi_annotation`) needs to be [included in your `$REFGENIE` configuration file](#Example_PEPPRO_REFGENIE_configuration_file) for the pipeline to detect it automatically.  Alternatively, you can use the `--pi-name` pipeline option to provide a path directly to this file.

### Premature mRNA

To determine the [*F*raction of *R*eads *i*n *P*re-mature mRNA (*FRiP*)](../glossary.md), you will need a [pre-mature mRNA annotation file](http://big.databio.org/peppro/). If a pre-built version for your genome of interest isn't present, you can create that file yourself. In the reference genome directory, execute the following (for `hg38`):
```console
wget -O hg38_refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
zcat hg38_refGene.txt.gz | grep 'cmpl' | \
  awk  '{print $3"\t"$5"\t"$6"\t"$4"\t"$13}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_pre-mRNA.tsv
```
This asset (`pre_mRNA_annotation`) needs to be [included in your `$REFGENIE` configuration file](#Example_PEPPRO_REFGENIE_configuration_file) for the pipeline to detect it automatically.  Alternatively, you can use the `--pre-name` pipeline option to provide a path directly to this file.

### Features

We also have [downloadable genome feature annotation files](http://big.databio.org/peppro/) for both `hg38` and `hg19` that you can use.  These files annotate 3' and 5' UTR, Exons, Introns, Promoters, and Promoter Flanking Regions.  If present in the corresponding reference genome folder and included as an asset (named `feat_annotation`) in your `$REFGENIE` configuration file you don't need to do anything else as the pipeline will look there automatically.   Alternatively, you can use the `--anno-name` pipeline option to just directly point to this file.

#### Create a custom feature annotation file

The pipeline will calculate the fraction of reads in genomic features using one of our [provided annotation files](http://big.databio.org/peppro/), but you can also specify this file yourself.

This annotation file is really just a modified `BED` file, with the chromosomal coordinates and type of feature included.  For example, the [downloadable `hg38_annotations.bed.gz` file](http://big.databio.org/peppro/hg38_annotations.bed.gz) looks like so:

```
chr1	28200	30001	Promoter	.	*
chr1	198800	200201	Promoter	.	*
chr1	778000	780001	Promoter	.	*
chr1	817400	817601	Promoter	.	*
chr1	826200	828801	Promoter	.	*
chr1	904200	905201	Promoter	.	*
chr1	923800	924601	Promoter	.	*
chr1	925000	925601	Promoter	.	*
chr1	941800	942201	Promoter	.	*
chr1	958400	961401	Promoter	.	*
```

Just like a standard `BED` file, the first three fields are:  
1. **chrom** - the name of the chromosome  
2. **chromStart** - the starting position of the feature  
3. **chromEnd** - the ending position of the feature

Column four is the **name** column, in our case the name of our feature of interest. The fifth column is the **score**, which would determine how darkly an item would be displayed in a genome browser if you chose to set that or if the information in your file of interest has ascribed a score to the features. The final, sixth, column is the **strand** column.

After creating your `BED` file, you can point the pipeline to it using the `--anno-name` option followed with the path to your file.  The pipeline will then use that file to determine the fractions of reads that cover those features.

### Example `PEPPRO` `REGENIE` configuration file

As mentioned above, you can point the pipeline directly to your annotation files using the matching arguments.

Alternatively, if they are all present in the corresponding reference genome folders, you can direct `refgenie` to detect them automatically. Here's an example of what a `refgenie` configuration file would look like:
```yaml
genome_folder: $GENOMES
genome_server: http://refgenomes.databio.org
genomes:
  hg38:
    bowtie2:
      path: indexed_bowtie2
    chrom_sizes:
      path: hg38.chrom.sizes
    tss_annotation:
      path: hg38_TSS.tsv
    cpa_annotation:
      path: hg38_CpA.tsv
    pre_mRNA_annotation:
      path: hg38_pre-mRNA.tsv
    feat_annotation:
      path: hg38_annotations.bed.gz
```