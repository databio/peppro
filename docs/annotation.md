# Annotation reference

This document outlines how we created the reference data, so you can recreate it if you need to. The easiest way to do this is use `refgenie build`. All you need to do is:


## 1: Build the fasta asset
You need a FASTA file for your genome. You can insert this file into refgenie like this:
```
refgenie build -g GENOME -a fasta --fasta path/to/file.fa
```

## 2: Build the GTF asset

You also need an Ensembl GTF file (or equivalent) for your genome. You can insert this file into refgenie like this:

```
refgenie build -g GENOME -a ensembl_gtf --GTF path/to/file.gtf
```

## 3: Build all other assets
Once you have those two assets installed, `refgenie` can automatically build all the remaining assets from them. Build the assets that are required like this:

```
refgenie build -g GENOME -a ensembl_gtf tss_annotation pre_mRNA_annotation \
  feat_annotation pi_tss pi_body exon_annotation intron_annotation
```

That's it! These assets will be automatically detected by PEPPRO if you build them like this with refgenie. If you want to know what we're doing, or customize these, more details follow:

### TSS

To calculate [TSS enrichments](../glossary.md), you will need a [TSS annotation file](http://big.databio.org/refgenomes/).  We build these using these commands, in this example for `hg38`:
```console
wget -O hg38_TSS_full.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz \
zcat hg38_TSS_full.txt.gz | \
  awk  '{if($4=="+"){print $3"\t"$5"\t"$5"\t"$13"\t.\t"$4}else{print $3"\t"$6"\t"$6"\t"$13"\t.\t"$4}}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_TSS.bed
```
You can pass the `--TSS-name` pipeline option to provide a path directly to this file.

### Pause index annotation (PI)

To calculate [pause indicies](../glossary.md), you will need two files in your reference genome directory: a PI TSS annotation file and a PI gene body annotation file. Here are example commands for `hg38`:
```console
wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.gtf.gz \
zcat Homo_sapiens.GRCh38.97.gtf.gz | \
  grep 'exon_number "1"' | \
  sed 's/^/chr/' | \
  awk '{OFS="\t";} {print $1,$4,$5,$20,$14,$7}' | \
  sed 's/";//g' | \
  sed 's/"//g' | \
  awk '{if($6=="+"){print $1"\t"$2+20"\t"$2+120"\t"$4"\t"$5"\t"$6}else{print $1"\t"$3-120"\t"$3-20"\t"$4"\t"$5"\t"$6}}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_PI_TSS.bed

zcat Homo_sapiens.GRCh38.97.gtf.gz | \
  awk '$3 == "gene"' | \
  sed 's/^/chr/' | \
  awk '{OFS="\t";} {print $1,$4,$5,$14,$6,$7}' | \
  sed 's/";//g' | \
  sed 's/"//g' |
  awk '$4!="Metazoa_SRP"' | \
  awk '$4!="U3"' | \
  awk '$4!="7SK"'  | \
  awk '($3-$2)>200' | \
  awk '{if($6=="+"){print $1"\t"$2+500"\t"$3"\t"$4"\t"$5"\t"$6}else{print $1"\t"$2"\t"$3-500"\t"$4"\t"$5"\t"$6}}' | \
  awk '$3>$2' | \
  LC_COLLATE=C sort -k4 -u > hg38_PI_gene_body.bed
```
You can use the `--pi-tss` and `--pi-body` pipeline options to provide paths directly to each file.

### mRNA contamination

To determine the amount of [mRNA contamination](../glossary.md), you will need two files in your reference genome directory: an [exon annotation file](http://big.databio.org/refgenomes/) and an [intron annotation file](http://big.databio.org/refgenomes/). Here are example commands for `hg38`:

```console
wget -O hg38_TSS_full.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz \
zcat hg38_TSS_full.txt.gz | \
  awk -v OFS="\t" '$9>1' | \
  awk -v OFS="\t" '{ n = split($10, a, ","); split($11, b, ","); for(i=1; i<n; ++i) print $3, a[i], b[i], $13, i, $4 }' | awk -v OFS="\t" '$6=="+" && $5!=1 {print $0} $6=="-" {print $0}' | \
  awk '$4!=prev4 && prev6=="-" {prev4=$4; prev6=$6; delete line[NR-1]; idx-=1} {line[++idx]=$0; prev4=$4; prev6=$6} END {for (x=1; x<=idx; x++) print line[x]}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_exons.bed

zcat hg38_TSS_full.txt.gz | \
  awk -v OFS="\t" '$9>1' | \
  awk -F"\t" '{ exonCount=int($9);split($10,exonStarts,"[,]"); split($11,exonEnds,"[,]"); for(i=1;i<exonCount;i++) {printf("%s\t%s\t%s\t%s\t%d\t%s\n",$3,exonEnds[i],exonStarts[i+1],$13,($3=="+"?i:exonCount-i),$4);}}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_introns.bed
```
You can use the `--exon-name` and `--intron-name` pipeline options to provide paths directly to each file.

### Premature mRNA

To determine the [*F*raction of *R*eads *i*n *P*re-mature mRNA (*FRiP*)](../glossary.md), you will need a [pre-mature mRNA annotation file](http://big.databio.org/peppro/). If a pre-built version for your genome of interest isn't present, you can create that file yourself. In the reference genome directory, execute the following (for `hg38`):
```console
wget -O hg38_refGene.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
zcat hg38_refGene.txt.gz | grep 'cmpl' | \
  awk  '{print $3"\t"$5"\t"$6"\t"$13"\t.\t"$4}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_pre-mRNA.bed
```
You can use the `--pre-name` pipeline option to provide a path directly to this file.

### Features

We also have [downloadable genome feature annotation files](http://big.databio.org/peppro/) for both `hg38` and `hg19` that you can use.  These files annotate 3' and 5' UTR, Exons, Introns, Promoters, and Promoter Flanking Regions.  You can use the `--anno-name` pipeline option to directly point to this file.

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
