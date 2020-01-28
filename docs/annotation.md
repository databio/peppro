# Custom reference data

The pipeline uses reference data at various stages, such as for alignment, calculating TSS enrichments, and other QC scores. If you're using a common genome assembly, these resources are pre-built and can be easily downloaded using `refgenie pull`, as described in the setup instructions. If the resources are not available, you'll have to build them. This document outlines how we created the reference data, so you can recreate it if you need to. The easiest way to do this is use `refgenie build`. All you need to do is:

## 1: Build the fasta asset

You need a FASTA file for your genome. You can insert this file into refgenie like this:
```console
refgenie build -g GENOME -a fasta --files fasta=/path/to/file.fa
```

## 2. Build the bowtie2_index

To build a bowtie2_index and have it managed by `refgenie` you'll, of course, need [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) already installed.  You will also need the requisite FASTA file, which you just added in step 1.
```console
refgenie build -g GENOME -a bowtie2_index
```

## 3: Build the ensembl_gtf asset

The ensembl_gtf asset includes several related assets (*e.g.* pause index gene bodies and TSS's) the pipeline will employ.  To build an ensembl_gtf asset, you need an Ensembl GTF file (or equivalent) for your genome. You can have refgenie build and manage this file as follows:

```console
refgenie build -g GENOME -a ensembl-gtf --files ensembl_gtf=/path/to/Homo_sapiens.GRCh38.97.gtf.gz
```

## 4: Build the refgene_anno asset

The refgene_anno asset actually includes several related assets that we'll need (*e.g.* TSS and premature mRNA annotations).  To build these, for example for hg38, you will need to [download a refGene annotation](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz). Build it for a any genome like so:

```console
refgenie build -g GENOME -a refgene_anno --files refgene=/path/to/refGene.txt.gz
```

## 5: Build the feat_annotation asset
The `feat_annotation` asset includes feature annotations used to calculate the [FRiF](glossary.md) and [PRiF](glossary.md). `Refgenie` can automatically build this after you have the above assets installed:

```console
refgenie build -g GENOME -a feat_annotation
```

That's it! These assets will be automatically detected by PEPPRO if you build them like this with `refgenie`. 

### Create a custom feature annotation file

The pipeline will calculate the fraction (and proportion) of reads in genomic features using the feat_annotation asset, but you can also construct this file yourself.

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
