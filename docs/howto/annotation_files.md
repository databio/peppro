# How to download or create annotation files

## Obtain annotation files

For each annotation type (TSS, premature mRNA, or general features), we provide [downloadable defaults](http://big.databio.org/peppro/) for common genomes.  You may also recreate these yourself as described below.

### TSS

To calculate TSS enrichments, you will need a [TSS annotation file](http://big.databio.org/refgenomes/) in your reference genome directory.  If a pre-built version for your genome of interest isn't present, you can quickly create that file yourself. In the reference genome directory, you can perform the following commands for in this example, `hg38`:
```console
wget -O hg38_TSS_full.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz \
zcat hg38_TSS_full.txt.gz | \
  awk  '{if($4=="+"){print $3"\t"$5"\t"$5"\t"$4"\t"$13}else{print $3"\t"$6"\t"$6"\t"$4"\t"$13}}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_TSS.tsv
```

### Premature mRNA

To determine the *F*raction of *R*eads *i*n *P*re-mature mRNA (*FRiP*), you will need a [pre-mature mRNA annotation file](http://big.databio.org/peppro/). If a pre-built version for your genome of interest isn't present, you can create that file yourself. In the reference genome directory, execute the following (for `hg38`):
```console
wget -O hg38_refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
zcat hg38_refGene.txt.gz | grep 'cmpl' | \
  awk  '{if($4=="+"){print $3"\t"$5"\t"$8"\t"$4"\t"$13}else{print $3"\t"$7"\t"$6"\t"$4"\t"$13}}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_pre-mRNA.tsv
```
The pipeline will look in the reference genome directory automatically.  Alternatively, you can use the `--pre-name` pipeline option to provide a path to this file.

### Features

We also have [downloadable genome feature annotation files](http://big.databio.org/peppro/) for both `hg38` and `hg19` that you can use.  These files annotate 3' and 5' UTR, Exons, Introns, Promoters, and Promoter Flanking Regions.  Simply move the corresponding genome annotation file into the `peppro/anno` folder.  Once present in the `peppro/anno` folder you don't need to do anything else as the pipeline will look there automatically.   Alternatively, you can use the `--anno-name` pipeline option to directly point to this file when running.

#### Create a custom feature annotation file

The pipeline will calculate the fraction of reads in genomic features using one of our [provided annotation files](http://big.databio.org/peppro/), but you can also specify this file yourself.

This annotation file is really just a modified `BED` file, with the chromosomal coordinates and type of feature included.  For example, the [downloadable `hg19_annotations.bed.gz` file](http://big.databio.org/peppro/hg19_annotations.bed.gz) looks like so:

```
chr1	28400	29801	Promoter	.	*
chr1	713400	715201	Promoter	.	*
chr1	761800	764001	Promoter	.	*
chr1	839600	840601	Promoter	.	*
chr1	859200	860001	Promoter	.	*
chr1	893800	895001	Promoter	.	*
chr1	895400	896801	Promoter	.	*
chr1	902000	902601	Promoter	.	*
chr1	934200	934801	Promoter	.	*
chr1	935200	936001	Promoter	.	*
```

Just like a standard `BED` file, the first three fields are:  
1. **chrom** - the name of the chromosome  
2. **chromStart** - the starting position of the feature  
3. **chromEnd** - the ending position of the feature

Column four is the **name** column, in our case the name of our feature of interest. The fifth column is the **score**, which would determine how darkly an item would be displayed in a genome browser if you chose to set that or if the information in your file of interest has ascribed a score to the features. The final, sixth, column is the **strand** column.

After creating your `BED` file, you can point the pipeline to it using the `--anno-name` option and providing the path to your file.  The pipeline will then use that file to determine the fractions of reads that cover those features.