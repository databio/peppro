### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_100 --genome hg38 --input /project/shefflab/data//sra_fastq/SRR1554311.fastq.gz /project/shefflab/data//sra_fastq/SRR1554312.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 32 -M 32000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aj40-15c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/
*  Pipeline started at:   (06-11 19:07:43) elapsed: 3.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.8
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16373 insertions(+), 3522 deletions(-)

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `cutadapt`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `32`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR1554311.fastq.gz', '/project/shefflab/data//sra_fastq/SRR1554312.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `32000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_PRO-seq_100`
*              `scale`:  `True`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `SINGLE`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `0`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /project/shefflab/data//sra_fastq/SRR1554311.fastq.gz

> `File_mb`	38578.5	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:07:44) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/K562_PRO-seq_100.merged.fastq.gz`  

> `cat /project/shefflab/data//sra_fastq/SRR1554311.fastq.gz /project/shefflab/data//sra_fastq/SRR1554312.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/K562_PRO-seq_100.merged.fastq.gz` (239583)
<pre>
</pre>
Command completed. Elapsed time: 0:02:48. Running peak memory: 0.003GB.  
  PID: 239583;	Command: cat;	Return code: 0;	Memory used: 0.003GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/K562_PRO-seq_100.merged.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1.fastq`  

> `pigz -f -p 32 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/K562_PRO-seq_100.merged.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1.fastq` (240015)
<pre>
</pre>
Command completed. Elapsed time: 0:14:07. Running peak memory: 0.003GB.  
  PID: 240015;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	496581677	PEPPRO	_RES_

> `Fastq_reads`	496581677	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/K562_PRO-seq_100.merged.fastq.gz']

### FASTQ processing:  (06-11 19:36:58) elapsed: 1755.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1_processed.fastq`  

> `(cutadapt -j 32 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt` (243484)
<pre>
</pre>
Command completed. Elapsed time: 0:14:01. Running peak memory: 11.493GB.  
  PID: 243484;	Command: cutadapt;	Return code: 0;	Memory used: 11.493GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1_processed.fastq` (245187,245188)
<pre>
</pre>
Command completed. Elapsed time: 0:09:01. Running peak memory: 11.493GB.  
  PID: 245187;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 245188;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	423059183.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	12437983.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.5047	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastqc/K562_PRO-seq_100_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (247019)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 11.493GB.  
  PID: 247019;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	484143694	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1_processed.fastq` (247474)
<pre>
Started analysis of K562_PRO-seq_100_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_100_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_100_R1_processed.fastq
Analysis complete for K562_PRO-seq_100_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:26:15. Running peak memory: 11.493GB.  
  PID: 247474;	Command: fastqc;	Return code: 0;	Memory used: 0.3GB

> `FastQC report r1`	fastqc/K562_PRO-seq_100_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/processed_R1.flag` (250743)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 11.493GB.  
  PID: 250743;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-11 20:34:31) elapsed: 3452.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt` (250749)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 11.493GB.  
  PID: 250749;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_100_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_100_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 20:34:37) elapsed: 6.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2312	PEPPRO	_RES_

### Prealignments (06-11 20:34:37) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 20:34:37) elapsed: 0.0 _TIME_


> `(bowtie2 -p 32 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq_100 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/prealignments/K562_PRO-seq_100_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
484143694 reads; of these:
  484143694 (100.00%) were unpaired; of these:
    439463330 (90.77%) aligned 0 times
    44680364 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	44680364.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-11 21:34:48) elapsed: 3611.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam`  

> `bowtie2 -p 32 --very-sensitive --rg-id K562_PRO-seq_100 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/prealignments/K562_PRO-seq_100_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/tmplun3n5k_ -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_temp.bam` (258021,258027,258028)
<pre>
439463330 reads; of these:
  439463330 (100.00%) were unpaired; of these:
    5393540 (1.23%) aligned 0 times
    328159166 (74.67%) aligned exactly 1 time
    105910624 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 142 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 3:21:41. Running peak memory: 11.493GB.  
  PID: 258021;	Command: bowtie2;	Return code: 0;	Memory used: 4.227GB  
  PID: 258027;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 258028;	Command: samtools;	Return code: 0;	Memory used: 0.906GB


> `samtools view -q 10 -b -@ 32 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam` (279311)
<pre>
</pre>
Command completed. Elapsed time: 0:08:44. Running peak memory: 11.493GB.  
  PID: 279311;	Command: samtools;	Return code: 0;	Memory used: 0.036GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	434069790	PEPPRO	_RES_

> `QC_filtered_reads`	46741044	PEPPRO	_RES_

> `Aligned_reads`	387328746	PEPPRO	_RES_

> `Alignment_rate`	80.0	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	29.61	PEPPRO	_RES_

### Compress all unmapped read files (06-12 02:01:41) elapsed: 16013.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/prealignments/K562_PRO-seq_100_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 32 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/prealignments/K562_PRO-seq_100_human_rDNA_unmap.fq` (286219)
<pre>
</pre>
Command completed. Elapsed time: 0:06:42. Running peak memory: 11.493GB.  
  PID: 286219;	Command: pigz;	Return code: 0;	Memory used: 0.024GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_temp.bam` (286889)
<pre>
</pre>
Command completed. Elapsed time: 0:06:33. Running peak memory: 11.493GB.  
  PID: 286889;	Command: samtools;	Return code: 0;	Memory used: 0.027GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	9149071	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam` (287540)
<pre>
</pre>
Command completed. Elapsed time: 0:05:38. Running peak memory: 11.493GB.  
  PID: 287540;	Command: samtools;	Return code: 0;	Memory used: 0.024GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/chr_sizes.bed` (288310,288311,288312,288313)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 11.493GB.  
  PID: 288312;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 288310;	Command: samtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 288313;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 288311;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/chr_sizes.bed -b -@ 32 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_noMT.bam` (288315)
<pre>
</pre>
Command completed. Elapsed time: 0:04:16. Running peak memory: 11.493GB.  
  PID: 288315;	Command: samtools;	Return code: 0;	Memory used: 0.037GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam` (288608)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 11.493GB.  
  PID: 288608;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam` (288612)
<pre>
</pre>
Command completed. Elapsed time: 0:05:31. Running peak memory: 11.493GB.  
  PID: 288612;	Command: samtools;	Return code: 0;	Memory used: 0.024GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-12 02:47:15) elapsed: 2733.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam -c 32 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_bamQC.tsv` (291051)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/tmp_K562_PRO-seq_100_sort_440xxpb4'
Processing with 32 cores...
Discarding 81 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1']
Keeping 114 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:04:01. Running peak memory: 16.596GB.  
  PID: 291051;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 16.596GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_bamQC.tsv`

> `NRF`	0.28	PEPPRO	_RES_

> `PBC1`	0.54	PEPPRO	_RES_

> `PBC2`	3.89	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_unmap.bam`  

> `samtools view -b -@ 32 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_unmap.bam` (291601)
<pre>
</pre>
Command completed. Elapsed time: 0:01:11. Running peak memory: 16.596GB.  
  PID: 291601;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `samtools view -c -f 4 -@ 32 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_temp.bam`

> `Unmapped_reads`	5393540	PEPPRO	_RES_

### Split BAM by strand (06-12 02:53:38) elapsed: 383.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam` (291818)
<pre>
</pre>
Command completed. Elapsed time: 0:27:15. Running peak memory: 16.596GB.  
  PID: 291818;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam` (294786)
<pre>
</pre>
Command completed. Elapsed time: 0:25:39. Running peak memory: 16.596GB.  
  PID: 294786;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-12 03:46:32) elapsed: 3174.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (297442)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 16.596GB.  
  PID: 297442;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/plus_TSS.tsv -p ends -c 32 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_plus_TssEnrichment.txt` (297448)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 16.596GB.  
  PID: 297448;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.744GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/minus_TSS.tsv -p ends -c 32 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_minus_TssEnrichment.txt` (297553)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 16.596GB.  
  PID: 297553;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 2.008GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_minus_TssEnrichment.txt` (297652)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 16.596GB.  
  PID: 297652;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_100_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_100_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt` (297684,297685,297686,297687)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 16.596GB.  
  PID: 297684;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 297686;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 297685;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 297687;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_keep.txt` (297689)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 16.596GB.  
  PID: 297689;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-12 03:47:43) elapsed: 71.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_ensembl_tss.bed` (297691,297692)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 16.596GB.  
  PID: 297691;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 297692;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_ensembl_gene_body.bed` (297696,297697)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 16.596GB.  
  PID: 297696;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 297697;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_TSS_density.bed` (297699,297700,297701,297702)
<pre>
</pre>
Command completed. Elapsed time: 0:08:08. Running peak memory: 16.596GB.  
  PID: 297699;	Command: bedtools;	Return code: 0;	Memory used: 0.103GB  
  PID: 297701;	Command: sort;	Return code: 0;	Memory used: 0.01GB  
  PID: 297700;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 297702;	Command: sort;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_gene_body_density.bed` (298620,298626,298627)
<pre>
</pre>
Command completed. Elapsed time: 0:14:05. Running peak memory: 16.596GB.  
  PID: 298626;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 298620;	Command: bedtools;	Return code: 0;	Memory used: 1.022GB  
  PID: 298627;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/tmpuv374kz8` (299969,299975,299976)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 16.596GB.  
  PID: 299969;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 299976;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 299975;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/tmpuv374kz8 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-12 04:09:58) elapsed: 1335.0 _TIME_

Total time: 9:02:18
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/tmpuv374kz8 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_100 --genome hg38 --input /project/shefflab/data//sra_fastq/SRR1554311.fastq.gz /project/shefflab/data//sra_fastq/SRR1554312.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 32 -M 32000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aj37-17c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/
*  Pipeline started at:   (06-14 21:10:54) elapsed: 6.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.8
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16373 insertions(+), 3522 deletions(-)

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `cutadapt`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `32`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR1554311.fastq.gz', '/project/shefflab/data//sra_fastq/SRR1554312.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `32000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_PRO-seq_100`
*              `scale`:  `True`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `SINGLE`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `0`
*          `verbosity`:  `None`

----------------------------------------

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//sra_fastq/SRR1554311.fastq.gz

> `File_mb`	38578.5	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:10:55) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/K562_PRO-seq_100.merged.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/K562_PRO-seq_100.merged.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/K562_PRO-seq_100_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/K562_PRO-seq_100.merged.fastq.gz']

### FASTQ processing:  (06-14 21:10:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:10:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/cutadapt/K562_PRO-seq_100_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_100_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_100_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:10:55) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:10:55) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/prealignments/K562_PRO-seq_100_human_rDNA_unmap.fq

### Map to genome (06-14 21:10:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam`  

### Compress all unmapped read files (06-14 21:10:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/prealignments/K562_PRO-seq_100_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:10:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_unmap.bam`  

### Split BAM by strand (06-14 21:10:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam`  

### Calculate TSS enrichment (06-14 21:10:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_PRO-seq_100_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_100_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:10:55) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/tmpot1hg31_` (280310,280311,280312)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 280310;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 280312;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 280311;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/tmpot1hg31_ | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.208447) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/tmpot1hg31_ > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_pause_index.bed` (280320)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 280320;	Command: awk;	Return code: 0;	Memory used: 0.004GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.96	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_pause_index.bed` (280325)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.212GB.  
  PID: 280325;	Command: Rscript;	Return code: 0;	Memory used: 0.212GB

> `Pause index`	QC_hg38/K562_PRO-seq_100_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_100_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_pause_index.bed.gz`  

> `pigz -f -p 32 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_pause_index.bed` (280358)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.212GB.  
  PID: 280358;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:11:01) elapsed: 6.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam`
387328746 138157932

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam`
387328746 131399421

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_gene_sort.bed` (331289,331339)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.212GB.  
  PID: 331289;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 331339;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_gene_coverage.bed` (331715)
<pre>
</pre>
Command completed. Elapsed time: 0:15:16. Running peak memory: 1.021GB.  
  PID: 331715;	Command: bedtools;	Return code: 0;	Memory used: 1.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/hg38_annotations.bed.gz` (393358)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.021GB.  
  PID: 393358;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 32 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/hg38_annotations.bed` (393364)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.021GB.  
  PID: 393364;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:37:32) elapsed: 1590.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/raw/hg38_annotations.bed` (393373)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 1.021GB.  
  PID: 393373;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Enhancer_sort.bed` (393375,393376,393377,393378)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 1.021GB.  
  PID: 393375;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 393376;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 393378;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB  
  PID: 393377;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Enhancer_plus_coverage.bed` (393381)
<pre>
</pre>
Command completed. Elapsed time: 0:04:20. Running peak memory: 1.021GB.  
  PID: 393381;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Enhancer_minus_coverage.bed` (393904)
<pre>
</pre>
Command completed. Elapsed time: 0:04:08. Running peak memory: 1.021GB.  
  PID: 393904;	Command: bedtools;	Return code: 0;	Memory used: 0.068GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_sort.bed` (394351,394352,394353,394354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.021GB.  
  PID: 394351;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 394352;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 394354;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 394353;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_plus_coverage.bed` (394356)
<pre>
</pre>
Command completed. Elapsed time: 0:05:01. Running peak memory: 1.021GB.  
  PID: 394356;	Command: bedtools;	Return code: 0;	Memory used: 0.592GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_minus_coverage.bed` (394913)
<pre>
</pre>
Command completed. Elapsed time: 0:04:57. Running peak memory: 1.021GB.  
  PID: 394913;	Command: bedtools;	Return code: 0;	Memory used: 0.475GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_Flanking_Region"` (395414)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.021GB.  
  PID: 395414;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_Flanking_Region_sort.bed` (395415,395416,395417,395418)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 1.021GB.  
  PID: 395415;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 395417;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 395416;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 395418;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_Flanking_Region_plus_coverage.bed` (395421)
<pre>
</pre>
Command completed. Elapsed time: 0:04:34. Running peak memory: 1.021GB.  
  PID: 395421;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_Flanking_Region_minus_coverage.bed` (395977)
<pre>
</pre>
Command completed. Elapsed time: 0:04:36. Running peak memory: 1.021GB.  
  PID: 395977;	Command: bedtools;	Return code: 0;	Memory used: 0.166GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/5_UTR"` (396492)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.021GB.  
  PID: 396492;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/5_UTR_sort.bed` (396493,396494,396495,396496)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 1.021GB.  
  PID: 396493;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 396494;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 396496;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 396495;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_5_UTR_plus_coverage.bed` (396500)
<pre>
</pre>
Command completed. Elapsed time: 0:04:30. Running peak memory: 1.021GB.  
  PID: 396500;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_5_UTR_minus_coverage.bed` (414548)
<pre>
</pre>
Command completed. Elapsed time: 0:04:12. Running peak memory: 1.021GB.  
  PID: 414548;	Command: bedtools;	Return code: 0;	Memory used: 0.063GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/3_UTR"` (419328)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.021GB.  
  PID: 419328;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/3_UTR_sort.bed` (419329,419330,419331,419332)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 1.021GB.  
  PID: 419329;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 419330;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 419332;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 419331;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_3_UTR_plus_coverage.bed` (419335)
<pre>
</pre>
Command completed. Elapsed time: 0:04:40. Running peak memory: 1.021GB.  
  PID: 419335;	Command: bedtools;	Return code: 0;	Memory used: 0.104GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_3_UTR_minus_coverage.bed` (419820)
<pre>
</pre>
Command completed. Elapsed time: 0:04:35. Running peak memory: 1.021GB.  
  PID: 419820;	Command: bedtools;	Return code: 0;	Memory used: 0.162GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Exon_sort.bed` (422032,422069,422085,422089)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 1.021GB.  
  PID: 422032;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 422069;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 422089;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 422085;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Exon_plus_coverage.bed` (424370)
<pre>
</pre>
Command completed. Elapsed time: 0:05:25. Running peak memory: 1.021GB.  
  PID: 424370;	Command: bedtools;	Return code: 0;	Memory used: 0.468GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Exon_minus_coverage.bed` (456033)
<pre>
</pre>
Command completed. Elapsed time: 0:05:15. Running peak memory: 1.021GB.  
  PID: 456033;	Command: bedtools;	Return code: 0;	Memory used: 0.186GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Intron_sort.bed` (29395,29396,29397,29398)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 1.021GB.  
  PID: 29395;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 29397;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 29396;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 29398;	Command: bedtools;	Return code: 0;	Memory used: 0.085GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Intron_plus_coverage.bed` (29402)
<pre>
</pre>
Command completed. Elapsed time: 0:06:15. Running peak memory: 1.021GB.  
  PID: 29402;	Command: bedtools;	Return code: 0;	Memory used: 0.453GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Intron_minus_coverage.bed` (57734)
<pre>
</pre>
Command completed. Elapsed time: 0:06:21. Running peak memory: 1.021GB.  
  PID: 57734;	Command: bedtools;	Return code: 0;	Memory used: 0.434GB


### Plot cFRiF/FRiF (06-14 22:46:31) elapsed: 4139.0 _TIME_


> `samtools view -@ 32 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_100 -z 3099922541 -n 192750289 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Intron_plus_coverage.bed` (80398)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 1.021GB.  
  PID: 80398;	Command: Rscript;	Return code: 0;	Memory used: 0.831GB

> `cFRiF`	QC_hg38/K562_PRO-seq_100_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_100_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_100 -z 3099922541 -n 192750289 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_Intron_plus_coverage.bed` (80459)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 1.021GB.  
  PID: 80459;	Command: Rscript;	Return code: 0;	Memory used: 0.439GB

> `FRiF`	QC_hg38/K562_PRO-seq_100_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_100_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 22:48:18) elapsed: 106.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_exons_sort.bed` (80502,80503)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 1.021GB.  
  PID: 80503;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB  
  PID: 80502;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_introns_sort.bed` (80512,80513,80514)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 1.021GB.  
  PID: 80512;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 80514;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 80513;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exons_coverage.bed` (80522)
<pre>
</pre>
Command completed. Elapsed time: 0:09:42. Running peak memory: 1.021GB.  
  PID: 80522;	Command: bedtools;	Return code: 0;	Memory used: 0.293GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_introns_coverage.bed` (131120)
<pre>
</pre>
Command completed. Elapsed time: 0:13:29. Running peak memory: 1.021GB.  
  PID: 131120;	Command: bedtools;	Return code: 0;	Memory used: 0.487GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/387.328746)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exons_rpkm.bed` (193263,193271,193272)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 1.021GB.  
  PID: 193263;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 193272;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 193271;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/387.328746)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_introns_rpkm.bed` (193274,193275,193276)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 1.021GB.  
  PID: 193274;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 193276;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 193275;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exon_intron_ratios.bed` (193279,193280,193281)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.021GB.  
  PID: 193279;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 193281;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 193280;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exon_intron_ratios.bed --annotate` (193287)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 1.021GB.  
  PID: 193287;	Command: Rscript;	Return code: 0;	Memory used: 0.208GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_100_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_100_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exon_intron_ratios.bed.gz`  

> `pigz -f -p 32 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/QC_hg38/K562_PRO-seq_100_exon_intron_ratios.bed` (193318)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.021GB.  
  PID: 193318;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (06-14 23:11:46) elapsed: 1408.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam` (193327)
<pre>
</pre>
Command completed. Elapsed time: 0:03:00. Running peak memory: 1.021GB.  
  PID: 193327;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_plus_smooth_body_0-mer.bw -p 21 --variable-step --tail-edge --scale 387328746.0` (193512)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_100_plus_cuttrace_l2uwmwym'
Processing with 10 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 90 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1']
Keeping 105 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 105 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_plus_exact_body_0-mer.bw'
Merging 105 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:30. Running peak memory: 3.847GB.  
  PID: 193512;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.847GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam` (195542)
<pre>
</pre>
Command completed. Elapsed time: 0:03:07. Running peak memory: 3.847GB.  
  PID: 195542;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_minus_smooth_body_0-mer.bw -p 21 --variable-step --tail-edge --scale 387328746.0` (195934)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/aligned_hg38/K562_PRO-seq_100_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_100_minus_cuttrace_l74fklbt'
Processing with 10 cores...
stdin is empty of data
psutil.NoSuchProcess process no longer exists (pid=195985)
Warning: couldn't add memory use for process: 195934
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
psutil.NoSuchProcess process no longer exists (pid=196723)
Warning: couldn't add memory use for process: 195934
Discarding 93 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 102 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270336v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 102 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_minus_exact_body_0-mer.bw'
Merging 102 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_100/signal_hg38/K562_PRO-seq_100_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:09. Running peak memory: 4.216GB.  
  PID: 195934;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 4.216GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:23:44
*  Total elapsed time (all runs):  13:29:45
*         Peak memory (this run):  4.2155 GB
*        Pipeline completed time: 2020-06-14 23:34:32
