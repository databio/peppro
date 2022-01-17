### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name Jurkat_ChRO-seq_1 --genome hg38 --input /project/shefflab/data//sra_fastq/SRR7616133.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 6 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-aj40-13c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/
*  Pipeline started at:   (06-11 17:35:31) elapsed: 1.0 _TIME_

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
*              `cores`:  `12`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR7616133.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `12000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `Jurkat_ChRO-seq_1`
*              `scale`:  `True`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `SINGLE`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `6`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /project/shefflab/data//sra_fastq/SRR7616133.fastq.gz

> `File_mb`	1330.25	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:35:33) elapsed: 2.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz`  

> `ln -sf /project/shefflab/data//sra_fastq/SRR7616133.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz` (26733)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 26733;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1.fastq` (26734)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 0.002GB.  
  PID: 26734;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	37740511	PEPPRO	_RES_

> `Fastq_reads`	37740511	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz']

### FASTQ processing:  (06-11 17:36:47) elapsed: 74.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_processed.fastq`  

> `(cutadapt -j 12 -m 8 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt` (26828)
<pre>
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 2.039GB.  
  PID: 26828;	Command: cutadapt;	Return code: 0;	Memory used: 2.039GB


> `seqtk trimfq -b 6 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_noadap.fastq | seqtk seq -L 8 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_processed.fastq` (26904,26905)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 2.039GB.  
  PID: 26904;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 26905;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	32783879	PEPPRO	_RES_

> `Trim_loss_rate`	13.13	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_processed.fastq` (26963)
<pre>
Started analysis of Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 5% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 10% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 15% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 20% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 25% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 30% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 35% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 40% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 45% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 50% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 55% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 60% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 65% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 70% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 75% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 80% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 85% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 90% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Approx 95% complete for Jurkat_ChRO-seq_1_R1_processed.fastq
Analysis complete for Jurkat_ChRO-seq_1_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 2.039GB.  
  PID: 26963;	Command: fastqc;	Return code: 0;	Memory used: 0.176GB

> `FastQC report r1`	fastqc/Jurkat_ChRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_noadap.fastq` (27048)
<pre>
[INFO][0m 14846147 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 2.039GB.  
  PID: 27048;	Command: seqkit;	Return code: 0;	Memory used: 1.036GB


> `seqtk trimfq -b 6 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_dedup.fastq | seqtk seq -L 8 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_trimmed.fastq` (27369,27370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 2.039GB.  
  PID: 27369;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 27370;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	34473620.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	2012355.0	PEPPRO	_RES_

> `Duplicate_reads`	14846147.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	5.3321	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/processed_R1.flag` (27478)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.039GB.  
  PID: 27478;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-11 17:41:41) elapsed: 293.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt -u 6` (27479)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 2.039GB.  
  PID: 27479;	Command: Rscript;	Return code: 0;	Memory used: 0.108GB

> `Adapter insertion distribution`	cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 17:41:52) elapsed: 12.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.9566	PEPPRO	_RES_

### Prealignments (06-11 17:41:53) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 17:41:53) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id Jurkat_ChRO-seq_1 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
32783879 reads; of these:
  32783879 (100.00%) were unpaired; of these:
    29290636 (89.34%) aligned 0 times
    3493243 (10.66%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.66% overall alignment rate

> `Aligned_reads_human_rDNA`	3493243.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.66	PEPPRO	_RES_

### Map to human_rDNA (06-11 17:44:54) elapsed: 181.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id Jurkat_ChRO-seq_1 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap_dups.fq 2>&1 > /dev/null)`

### Map to genome (06-11 17:46:43) elapsed: 109.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id Jurkat_ChRO-seq_1 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/tmpgdqhj09o -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam` (28061,28062,28063)
<pre>
29290636 reads; of these:
  29290636 (100.00%) were unpaired; of these:
    1113286 (3.80%) aligned 0 times
    17309978 (59.10%) aligned exactly 1 time
    10867372 (37.10%) aligned >1 times
96.20% overall alignment rate
[bam_sort_core] merging from 8 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:11:07. Running peak memory: 3.685GB.  
  PID: 28061;	Command: bowtie2;	Return code: 0;	Memory used: 3.685GB  
  PID: 28062;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 28063;	Command: samtools;	Return code: 0;	Memory used: 0.9GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam` (3548)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 3.685GB.  
  PID: 3548;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	28177350	PEPPRO	_RES_

> `QC_filtered_reads`	6839271	PEPPRO	_RES_

> `Aligned_reads`	21338079	PEPPRO	_RES_

> `Alignment_rate`	65.09	PEPPRO	_RES_

> `Total_efficiency`	56.54	PEPPRO	_RES_

> `Read_depth`	3.87	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id Jurkat_ChRO-seq_1 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap_dups.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/tmpgdqhj09o -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp_dups.bam` (4698,4703,4704)
<pre>
18024676 reads; of these:
  18024676 (100.00%) were unpaired; of these:
    856458 (4.75%) aligned 0 times
    10720803 (59.48%) aligned exactly 1 time
    6447415 (35.77%) aligned >1 times
95.25% overall alignment rate
[bam_sort_core] merging from 5 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:07:33. Running peak memory: 3.685GB.  
  PID: 4698;	Command: bowtie2;	Return code: 0;	Memory used: 3.684GB  
  PID: 4703;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 4704;	Command: samtools;	Return code: 0;	Memory used: 0.899GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam` (5505)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.685GB.  
  PID: 5505;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-11 18:14:58) elapsed: 1695.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq` (5562)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 3.685GB.  
  PID: 5562;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam` (5831)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.685GB.  
  PID: 5831;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	4738	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam` (5859)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.685GB.  
  PID: 5859;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/chr_sizes.bed` (5883,5884,5885,5886)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.685GB.  
  PID: 5883;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 5885;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 5884;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 5886;	Command: grep;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_noMT.bam` (5888)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.685GB.  
  PID: 5888;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam` (5926)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.685GB.  
  PID: 5926;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam` (5927)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.685GB.  
  PID: 5927;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	70	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp_dups.bam` (6012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.685GB.  
  PID: 6012;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam.bai

### Calculate library complexity (06-11 18:17:51) elapsed: 173.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam` (6027)
<pre>
BAM_INPUT
TOTAL READS     = 13139331
COUNTS_SUM      = 13139331
DISTINCT READS  = 9.95264e+06
DISTINCT COUNTS = 149
MAX COUNT       = 430
COUNTS OF 1     = 8.33143e+06
OBSERVED COUNTS (431)
1	8331431
2	1064073
3	273280
4	110730
5	57100
6	33562
7	21886
8	14351
9	10164
10	7277
11	5449
12	4060
13	3259
14	2572
15	2117
16	1690
17	1348
18	1122
19	924
20	809
21	697
22	515
23	473
24	408
25	362
26	314
27	268
28	197
29	186
30	150
31	147
32	163
33	117
34	133
35	90
36	92
37	78
38	88
39	59
40	50
41	60
42	49
43	34
44	46
45	41
46	32
47	44
48	36
49	27
50	22
51	33
52	22
53	21
54	16
55	20
56	17
57	20
58	12
59	14
60	10
61	9
62	8
63	12
64	9
65	11
66	6
67	8
68	5
69	10
70	8
71	8
72	7
73	3
74	8
75	3
76	6
77	5
78	2
79	3
80	6
81	5
82	6
83	7
84	3
85	4
86	2
87	2
88	3
89	3
90	4
91	2
92	6
93	2
94	4
95	3
96	2
98	1
99	2
100	2
101	1
102	3
104	1
105	1
107	2
108	2
109	1
110	5
111	2
112	1
113	1
114	1
116	2
117	2
118	2
123	1
125	1
126	1
128	1
132	1
133	2
134	1
135	1
136	1
139	1
141	1
143	1
144	1
147	2
150	1
151	1
156	2
167	2
168	1
172	2
174	1
177	1
189	1
191	1
193	2
215	1
216	1
225	1
280	1
296	1
349	1
361	1
388	1
402	1
430	1

sample size: 1000000
sample size: 2000000
sample size: 3000000
sample size: 4000000
sample size: 5000000
sample size: 6000000
sample size: 7000000
sample size: 8000000
sample size: 9000000
sample size: 10000000
sample size: 11000000
sample size: 12000000
sample size: 13000000
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 3.685GB.  
  PID: 6027;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam` (6135)
<pre>
BAM_INPUT
TOTAL READS     = 13139331
DISTINCT READS  = 9.95264e+06
DISTINCT COUNTS = 149
MAX COUNT       = 430
COUNTS OF 1     = 8.33143e+06
MAX TERMS       = 96
OBSERVED COUNTS (431)
1	8331431
2	1064073
3	273280
4	110730
5	57100
6	33562
7	21886
8	14351
9	10164
10	7277
11	5449
12	4060
13	3259
14	2572
15	2117
16	1690
17	1348
18	1122
19	924
20	809
21	697
22	515
23	473
24	408
25	362
26	314
27	268
28	197
29	186
30	150
31	147
32	163
33	117
34	133
35	90
36	92
37	78
38	88
39	59
40	50
41	60
42	49
43	34
44	46
45	41
46	32
47	44
48	36
49	27
50	22
51	33
52	22
53	21
54	16
55	20
56	17
57	20
58	12
59	14
60	10
61	9
62	8
63	12
64	9
65	11
66	6
67	8
68	5
69	10
70	8
71	8
72	7
73	3
74	8
75	3
76	6
77	5
78	2
79	3
80	6
81	5
82	6
83	7
84	3
85	4
86	2
87	2
88	3
89	3
90	4
91	2
92	6
93	2
94	4
95	3
96	2
98	1
99	2
100	2
101	1
102	3
104	1
105	1
107	2
108	2
109	1
110	5
111	2
112	1
113	1
114	1
116	2
117	2
118	2
123	1
125	1
126	1
128	1
132	1
133	2
134	1
135	1
136	1
139	1
141	1
143	1
144	1
147	2
150	1
151	1
156	2
167	2
168	1
172	2
174	1
177	1
189	1
191	1
193	2
215	1
216	1
225	1
280	1
296	1
349	1
361	1
388	1
402	1
430	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
.................._...._..............................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:13. Running peak memory: 3.685GB.  
  PID: 6135;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort_dups.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_counts.txt` (6527)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.685GB.  
  PID: 6527;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_plot` (6549)
<pre>
Processing Jurkat_ChRO-seq_1
INFO: Found real counts for Jurkat_ChRO-seq_1 - Total (M): 21.334642 Unique (M): 13.139331

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.685GB.  
  PID: 6549;	Command: Rscript;	Return code: 0;	Memory used: 0.085GB

> `Library complexity`	QC_hg38/Jurkat_ChRO-seq_1_preseq_plot.pdf	Library complexity	QC_hg38/Jurkat_ChRO-seq_1_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.7897	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 18:20:51) elapsed: 180.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv` (6575)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/tmp_Jurkat_ChRO-seq_1_sort_qn69eg_r'
Processing with 12 cores...
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.685GB.  
  PID: 6575;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.521GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv`

> `NRF`	0.38	PEPPRO	_RES_

> `PBC1`	0.64	PEPPRO	_RES_

> `PBC2`	3.39	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_unmap.bam` (6641)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.685GB.  
  PID: 6641;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_temp.bam`

> `Unmapped_reads`	1113286	PEPPRO	_RES_

### Split BAM by strand (06-11 18:21:24) elapsed: 33.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam` (6680)
<pre>
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 3.685GB.  
  PID: 6680;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam` (6747)
<pre>
</pre>
Command completed. Elapsed time: 0:00:59. Running peak memory: 3.685GB.  
  PID: 6747;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 18:23:23) elapsed: 119.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (6824)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.685GB.  
  PID: 6824;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_plus_TssEnrichment.txt` (6825)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.685GB.  
  PID: 6825;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.534GB


> `TSS_coding_score`	55.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_minus_TssEnrichment.txt` (6860)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.685GB.  
  PID: 6860;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.726GB


> `TSS_non-coding_score`	12.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_minus_TssEnrichment.txt` (6893)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.685GB.  
  PID: 6893;	Command: Rscript;	Return code: 0;	Memory used: 0.108GB

> `TSS enrichment`	QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.pdf	TSS enrichment	QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt` (6919,6920,6921,6922)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.685GB.  
  PID: 6919;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 6921;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 6920;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 6922;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt` (6924)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.685GB.  
  PID: 6924;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-11 18:23:49) elapsed: 26.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_tss.bed` (6926,6927)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.685GB.  
  PID: 6926;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 6927;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed` (6931,6932)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.685GB.  
  PID: 6931;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 6932;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSS_density.bed` (6934,6935,6936,6937)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.685GB.  
  PID: 6935;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 6937;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 6934;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB  
  PID: 6936;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_gene_body_density.bed` (6977,6978,6979)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 3.685GB.  
  PID: 6979;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 6977;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB  
  PID: 6978;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed`  
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:             ^ syntax error
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:                                                                                                                      ^ syntax error

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/tmpf89gr92g` (7009,7010,7011)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:14) elapsed: 1885.0 _TIME_
Child process 7009 (join) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/recover.lock.QC_hg38__Jurkat_ChRO-seq_1_pause_index.bed

### Pipeline failed at:  (06-11 18:55:14) elapsed: 0.0 _TIME_

Total time: 1:19:45
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name Jurkat_ChRO-seq_1 --genome hg38 --input /project/shefflab/data//sra_fastq/SRR7616133.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 6 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba26-12
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/
*  Pipeline started at:   (06-11 19:11:35) elapsed: 0.0 _TIME_

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
*              `cores`:  `12`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR7616133.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `12000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `Jurkat_ChRO-seq_1`
*              `scale`:  `True`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `SINGLE`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `6`
*          `verbosity`:  `None`

----------------------------------------

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//sra_fastq/SRR7616133.fastq.gz

> `File_mb`	1330.25	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:11:36) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz']

### FASTQ processing:  (06-11 19:11:36) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-11 19:11:36) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-11 19:11:36) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:11:36) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq

### Map to human_rDNA (06-11 19:11:36) elapsed: 0.0 _TIME_


### Map to genome (06-11 19:11:36) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam`  

### Compress all unmapped read files (06-11 19:11:36) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-11 19:11:36) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_unmap.bam`  

### Split BAM by strand (06-11 19:11:36) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam`  

### Calculate TSS enrichment (06-11 19:11:36) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.pdf	TSS enrichment	QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-11 19:11:36) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_gene_body_density.bed`  
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/lock.QC_hg38__Jurkat_ChRO-seq_1_pause_index.bed
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/tmpztt2wubq` (45321,45322,45323)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 45321;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 45323;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 45322;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/tmpztt2wubq | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:11:36) elapsed: 0.0 _TIME_

Total time: 0:00:01
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/tmpztt2wubq | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name Jurkat_ChRO-seq_1 --genome hg38 --input /project/shefflab/data//sra_fastq/SRR7616133.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 6 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-28c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/
*  Pipeline started at:   (06-14 21:14:42) elapsed: 6.0 _TIME_

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
*              `cores`:  `12`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//sra_fastq/SRR7616133.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `12000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `Jurkat_ChRO-seq_1`
*              `scale`:  `True`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `SINGLE`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `6`
*          `verbosity`:  `None`

----------------------------------------

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//sra_fastq/SRR7616133.fastq.gz

> `File_mb`	1330.25	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:14:43) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/Jurkat_ChRO-seq_1_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/Jurkat_ChRO-seq_1.fastq.gz']

### FASTQ processing:  (06-14 21:14:43) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:14:43) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/Jurkat_ChRO-seq_1_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:14:43) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:14:43) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq

### Map to human_rDNA (06-14 21:14:43) elapsed: 0.0 _TIME_


### Map to genome (06-14 21:14:43) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam`  

### Compress all unmapped read files (06-14 21:14:43) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/prealignments/Jurkat_ChRO-seq_1_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:14:43) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_unmap.bam`  

### Split BAM by strand (06-14 21:14:43) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam`  

### Calculate TSS enrichment (06-14 21:14:43) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.pdf	TSS enrichment	QC_hg38/Jurkat_ChRO-seq_1_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:14:43) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/tmp6gr4kvrx` (304870,304873,304875)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 304870;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 304875;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 304873;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/tmp6gr4kvrx | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0149084) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/tmp6gr4kvrx > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed` (304905)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 304905;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	42.69	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed` (304910)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.244GB.  
  PID: 304910;	Command: Rscript;	Return code: 0;	Memory used: 0.244GB

> `Pause index`	QC_hg38/Jurkat_ChRO-seq_1_pause_index.pdf	Pause index	QC_hg38/Jurkat_ChRO-seq_1_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_pause_index.bed` (304933)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.244GB.  
  PID: 304933;	Command: pigz;	Return code: 0;	Memory used: 0.004GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:14:49) elapsed: 7.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam`
21338079 7688297

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam`
21338079 7358614

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_gene_sort.bed` (305367,305368)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.244GB.  
  PID: 305367;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 305368;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_gene_coverage.bed` (305370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 0.244GB.  
  PID: 305370;	Command: bedtools;	Return code: 0;	Memory used: 0.085GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed.gz` (305399)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.244GB.  
  PID: 305399;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed` (305400)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.244GB.  
  PID: 305400;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:15:57) elapsed: 68.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/raw/hg38_annotations.bed` (305409)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.244GB.  
  PID: 305409;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer_sort.bed` (305412,305413,305414,305415)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.244GB.  
  PID: 305412;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 305413;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 305415;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB  
  PID: 305414;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_plus_coverage.bed` (305417)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 0.244GB.  
  PID: 305417;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_minus_coverage.bed` (305428)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 0.244GB.  
  PID: 305428;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_sort.bed` (305538,305539,305540,305541)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.244GB.  
  PID: 305538;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 305539;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 305541;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 305540;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_plus_coverage.bed` (305554)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.244GB.  
  PID: 305554;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_minus_coverage.bed` (305661)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.244GB.  
  PID: 305661;	Command: bedtools;	Return code: 0;	Memory used: 0.031GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region"` (305675)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.244GB.  
  PID: 305675;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed` (305676,305677,305678,305679)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.244GB.  
  PID: 305676;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 305678;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 305677;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 305679;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (305682)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.244GB.  
  PID: 305682;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed` (307252)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.244GB.  
  PID: 307252;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR"` (316865)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.244GB.  
  PID: 316865;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR_sort.bed` (317008,317009,317010,317011)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.244GB.  
  PID: 317008;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 317009;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 317011;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 317010;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_plus_coverage.bed` (317625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.244GB.  
  PID: 317625;	Command: bedtools;	Return code: 0;	Memory used: 0.028GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_minus_coverage.bed` (332877)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.244GB.  
  PID: 332877;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR"` (348684)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.244GB.  
  PID: 348684;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR_sort.bed` (348882,348918,348937,348945)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.244GB.  
  PID: 348882;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 348918;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 348945;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 348937;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_plus_coverage.bed` (350591)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.244GB.  
  PID: 350591;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_minus_coverage.bed` (374792)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 0.244GB.  
  PID: 374792;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon_sort.bed` (397491,397549,397556,397558)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.244GB.  
  PID: 397491;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 397549;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 397558;	Command: bedtools;	Return code: 0;	Memory used: 0.177GB  
  PID: 397556;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_plus_coverage.bed` (403618)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 0.244GB.  
  PID: 403618;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_minus_coverage.bed` (429156)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 0.244GB.  
  PID: 429156;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron_sort.bed` (447116,447165,447182,447189)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.244GB.  
  PID: 447116;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 447182;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 447165;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 447189;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_plus_coverage.bed` (449088)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 0.244GB.  
  PID: 449088;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_minus_coverage.bed` (1087)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 0.244GB.  
  PID: 1087;	Command: bedtools;	Return code: 0;	Memory used: 0.049GB


### Plot cFRiF/FRiF (06-14 21:19:07) elapsed: 190.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s Jurkat_ChRO-seq_1 -z 3099922541 -n 10780321 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_plus_coverage.bed` (1141)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 0.469GB.  
  PID: 1141;	Command: Rscript;	Return code: 0;	Memory used: 0.469GB

> `cFRiF`	QC_hg38/Jurkat_ChRO-seq_1_cFRiF.pdf	cFRiF	QC_hg38/Jurkat_ChRO-seq_1_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s Jurkat_ChRO-seq_1 -z 3099922541 -n 10780321 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_Intron_plus_coverage.bed` (1220)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 0.469GB.  
  PID: 1220;	Command: Rscript;	Return code: 0;	Memory used: 0.467GB

> `FRiF`	QC_hg38/Jurkat_ChRO-seq_1_FRiF.pdf	FRiF	QC_hg38/Jurkat_ChRO-seq_1_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:20:10) elapsed: 63.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_exons_sort.bed` (1653,1654)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.469GB.  
  PID: 1654;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB  
  PID: 1653;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_introns_sort.bed` (1668,1669,1670)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.469GB.  
  PID: 1668;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 1670;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 1669;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_coverage.bed` (1694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 0.469GB.  
  PID: 1694;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_coverage.bed` (1756)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 0.469GB.  
  PID: 1756;	Command: bedtools;	Return code: 0;	Memory used: 0.068GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/21.338079)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_rpkm.bed` (1778,1779,1780)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.469GB.  
  PID: 1778;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 1780;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 1779;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/21.338079)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_rpkm.bed` (1782,1783,1784)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.469GB.  
  PID: 1782;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 1784;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 1783;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed` (1786,1787,1788)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.469GB.  
  PID: 1786;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 1788;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 1787;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.19	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed --annotate` (1795)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.469GB.  
  PID: 1795;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/Jurkat_ChRO-seq_1_mRNA_contamination.pdf	mRNA contamination	QC_hg38/Jurkat_ChRO-seq_1_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/QC_hg38/Jurkat_ChRO-seq_1_exon_intron_ratios.bed` (1815)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.469GB.  
  PID: 1815;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-14 21:21:13) elapsed: 63.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam` (1823)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.469GB.  
  PID: 1823;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 21338079.0` (1830)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_plus.bam'
Temporary files will be stored in: 'tmp_Jurkat_ChRO-seq_1_plus_cuttrace_gbzfszua'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 122 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr15_KI270727v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_KI270743v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270742v1']
Keeping 73 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270579v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_GL000213v1', 'chrUn_KI270744v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 73 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_exact_body_0-mer.bw'
Merging 73 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:37. Running peak memory: 2.442GB.  
  PID: 1830;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.442GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam` (311235)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 2.442GB.  
  PID: 311235;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 21338079.0` (311251)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/aligned_hg38/Jurkat_ChRO-seq_1_minus.bam'
Temporary files will be stored in: 'tmp_Jurkat_ChRO-seq_1_minus_cuttrace_or0c10qu'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 126 chunk(s) of reads: ['chrM', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr14_KI270725v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1']
Keeping 69 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270580v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 69 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_minus_exact_body_0-mer.bw'
Merging 69 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/Jurkat_ChRO-seq_1/signal_hg38/Jurkat_ChRO-seq_1_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:27. Running peak memory: 3.351GB.  
  PID: 311251;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.351GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:19:54
*  Total elapsed time (all runs):  1:16:55
*         Peak memory (this run):  3.3511 GB
*        Pipeline completed time: 2020-06-14 21:34:30
