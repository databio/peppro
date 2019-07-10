# Download or create annotation files for <img src="../../img/peppro_logo.svg" alt="PEPPRO" class="img-fluid" style="max-height:35px; margin-top:-15px; margin-bottom:-10px"> 


For each annotation type (TSS, CpA sites, premature mRNA, or general features), we provide [downloadable defaults](http://big.databio.org/peppro/) for common genomes.  You may also recreate these yourself as described below.

### TSS

To calculate [TSS enrichments](../glossary.md), you will need a [TSS annotation file](http://big.databio.org/refgenomes/) in your reference genome directory.  If a pre-built version for your genome of interest isn't present, you can quickly create that file yourself. In the reference genome directory, you can perform the following commands for in this example, `hg38`:
```console
wget -O hg38_TSS_full.txt.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz \
zcat hg38_TSS_full.txt.gz | \
  awk  '{if($4=="+"){print $3"\t"$5"\t"$5"\t"$4"\t"$13}else{print $3"\t"$6"\t"$6"\t"$4"\t"$13}}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_TSS.tsv
```
This asset (`tss_annotation`) needs to be [included in your `$REFGENIE` configuration file](annotation_files.md#example-peppro-refgenie-configuration-file) for the pipeline to detect it automatically.  Alternatively, you can use the `--TSS-name` pipeline option to provide a path directly to this file.

### Pause index annotation (PI)

To calculate [pause indicies](../glossary.md), you will need two files in your reference genome directory: a [PI TSS annotation file](http://big.databio.org/refgenomes/) and a [PI gene body annotation file](http://big.databio.org/refgenomes/).  If a pre-built version for your genome of interest isn't present, you can quickly create that file yourself. In the reference genome directory, you can perform the following commands for in this example, `hg38`:
```console
wget ftp://ftp.ensembl.org/pub/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh38.87.gtf.gz \
zcat Homo_sapiens.GRCh38.87.gtf.gz | \
  grep 'exon_number "1"' | \
  sed 's/^/chr/' | \
  awk '{OFS="\t";} {print $1,$4,$5,$20,$14,$7}' | \
  sed 's/";//g' | \
  sed 's/"//g' | \
  awk '{if($6=="+"){print $1"\t"$2+20"\t"$3+120"\t"$4"\t"$5"\t"$6}else{print $1"\t"$3-120"\t"$3-20"\t"$4"\t"$5"\t"$6}}' | \
  LC_COLLATE=C sort -k1,1 -k2,2n -u > hg38_PI_TSS.bed

zcat Homo_sapiens.GRCh38.87.gtf.gz | \
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
These assets (`pi_tss` and `pi_body`) need to be [included in your `$REFGENIE` configuration file](annotation_files.md#example-peppro-refgenie-configuration-file) for the pipeline to detect it automatically.  Alternatively, you can use the `--pi-tss` and `--pi-body` pipeline options to provide paths directly to each file.

### mRNA contamination

To determine the amount of [mRNA contamination](../glossary.md), you will need two files in your reference genome directory: an [exon annotation file](http://big.databio.org/refgenomes/) and an [intron annotation file](http://big.databio.org/refgenomes/).  If a pre-built version for your genome of interest isn't present, you can quickly create that file yourself. In the reference genome directory, you can perform the following commands for in this example, `hg38`:
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
These assets (`exon_annotation` and `intron_annotation`) need to be [included in your `$REFGENIE` configuration file](annotation_files.md#example-peppro-refgenie-configuration-file) for the pipeline to detect it automatically.  Alternatively, you can use the `--exon-name` and `--intron-name` pipeline options to provide paths directly to each file.

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

### Example `PEPPRO` `refgenie` configuration file

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
    pi_tss:
      path: hg38_PI_TSS.bed
    pi_body:
      path: hg38_PI_gene_body.bed
    pre_mRNA_annotation:
      path: hg38_pre-mRNA.tsv
    feat_annotation:
      path: hg38_annotations.bed.gz
    exon_annotation:
      path: hg38_exons.bed
    intron_annotation:
      path: hg38_introns.bed
```