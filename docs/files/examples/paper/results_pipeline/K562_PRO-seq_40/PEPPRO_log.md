### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_40 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_40pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 16 -M 16000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-26b
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/
*  Pipeline started at:   (06-15 07:11:05) elapsed: 4.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.8
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16413 insertions(+), 3702 deletions(-)

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `cutadapt`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `16`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_40pct.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `16000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_PRO-seq_40`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_40pct.fastq.gz

> `File_mb`	14593.59	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:11:05) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/K562_PRO-seq_40.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_PRO_40pct.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/K562_PRO-seq_40.fastq.gz` (208283)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 208283;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/K562_PRO-seq_40.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1.fastq`  

> `pigz -f -p 16 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/K562_PRO-seq_40.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1.fastq` (208285)
<pre>
</pre>
Command completed. Elapsed time: 0:08:26. Running peak memory: 0.003GB.  
  PID: 208285;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	198648534	PEPPRO	_RES_

> `Fastq_reads`	198648534	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/K562_PRO-seq_40.fastq.gz']

### FASTQ processing:  (06-15 07:28:19) elapsed: 1034.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1_processed.fastq`  

> `(cutadapt -j 16 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt` (209838)
<pre>
</pre>
Command completed. Elapsed time: 0:09:00. Running peak memory: 6.657GB.  
  PID: 209838;	Command: cutadapt;	Return code: 0;	Memory used: 6.657GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1_processed.fastq` (210820,210821)
<pre>
</pre>
Command completed. Elapsed time: 0:06:11. Running peak memory: 6.657GB.  
  PID: 210821;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB  
  PID: 210820;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	169233215.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	4973679.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.5038	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastqc/K562_PRO-seq_40_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (211947)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 211947;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	193674855	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1_processed.fastq` (212279)
<pre>
Started analysis of K562_PRO-seq_40_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_40_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_40_R1_processed.fastq
Analysis complete for K562_PRO-seq_40_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:11:03. Running peak memory: 6.657GB.  
  PID: 212279;	Command: fastqc;	Return code: 0;	Memory used: 0.241GB

> `FastQC report r1`	fastqc/K562_PRO-seq_40_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/processed_R1.flag` (213350)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 213350;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-15 08:02:23) elapsed: 2044.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt` (213355)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.657GB.  
  PID: 213355;	Command: Rscript;	Return code: 0;	Memory used: 0.239GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_40_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_40_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 08:02:29) elapsed: 6.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/cutadapt/K562_PRO-seq_40_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2312	PEPPRO	_RES_

### Prealignments (06-15 08:02:29) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 08:02:29) elapsed: 0.0 _TIME_


> `(bowtie2 -p 16 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq_40 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/fastq/K562_PRO-seq_40_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/prealignments/K562_PRO-seq_40_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
193674855 reads; of these:
  193674855 (100.00%) were unpaired; of these:
    175801164 (90.77%) aligned 0 times
    17873691 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	17873691.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-15 08:24:20) elapsed: 1312.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam`  

> `bowtie2 -p 16 --very-sensitive --rg-id K562_PRO-seq_40 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/prealignments/K562_PRO-seq_40_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/tmpk69t81xu -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_temp.bam` (216329,216334,216335)
<pre>
175801164 reads; of these:
  175801164 (100.00%) were unpaired; of these:
    2158413 (1.23%) aligned 0 times
    131278595 (74.67%) aligned exactly 1 time
    42364156 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 56 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 1:39:26. Running peak memory: 6.657GB.  
  PID: 216329;	Command: bowtie2;	Return code: 0;	Memory used: 3.852GB  
  PID: 216334;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 216335;	Command: samtools;	Return code: 0;	Memory used: 0.96GB


> `samtools view -q 10 -b -@ 16 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam` (231982)
<pre>
</pre>
Command completed. Elapsed time: 0:04:44. Running peak memory: 6.657GB.  
  PID: 231982;	Command: samtools;	Return code: 0;	Memory used: 0.022GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	173642751	PEPPRO	_RES_

> `QC_filtered_reads`	18693246	PEPPRO	_RES_

> `Aligned_reads`	154949505	PEPPRO	_RES_

> `Alignment_rate`	80.0	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	13.26	PEPPRO	_RES_

### Compress all unmapped read files (06-15 10:44:54) elapsed: 8434.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/prealignments/K562_PRO-seq_40_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/prealignments/K562_PRO-seq_40_human_rDNA_unmap.fq` (236098)
<pre>
</pre>
Command completed. Elapsed time: 0:04:32. Running peak memory: 6.657GB.  
  PID: 236098;	Command: pigz;	Return code: 0;	Memory used: 0.014GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_temp.bam` (236559)
<pre>
</pre>
Command completed. Elapsed time: 0:02:47. Running peak memory: 6.657GB.  
  PID: 236559;	Command: samtools;	Return code: 0;	Memory used: 0.024GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	3660275	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam` (236990)
<pre>
</pre>
Command completed. Elapsed time: 0:02:24. Running peak memory: 6.657GB.  
  PID: 236990;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/chr_sizes.bed` (237141,237142,237143,237144)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 237143;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 237141;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 237144;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 237142;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/chr_sizes.bed -b -@ 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_noMT.bam` (237146)
<pre>
</pre>
Command completed. Elapsed time: 0:02:45. Running peak memory: 6.657GB.  
  PID: 237146;	Command: samtools;	Return code: 0;	Memory used: 0.023GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam` (237517)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.657GB.  
  PID: 237517;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam` (237519)
<pre>
</pre>
Command completed. Elapsed time: 0:02:26. Running peak memory: 6.657GB.  
  PID: 237519;	Command: samtools;	Return code: 0;	Memory used: 0.017GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 11:07:25) elapsed: 1351.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam -c 16 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_bamQC.tsv` (238553)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/tmp_K562_PRO-seq_40_sort_j78vxxnp'
Processing with 16 cores...
Discarding 87 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 108 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:02:34. Running peak memory: 6.657GB.  
  PID: 238553;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 4.815GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_bamQC.tsv`

> `NRF`	0.53	PEPPRO	_RES_

> `PBC1`	0.74	PEPPRO	_RES_

> `PBC2`	5.58	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_unmap.bam`  

> `samtools view -b -@ 16 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_unmap.bam` (238762)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 6.657GB.  
  PID: 238762;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools view -c -f 4 -@ 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_temp.bam`

> `Unmapped_reads`	2158413	PEPPRO	_RES_

### Split BAM by strand (06-15 11:10:57) elapsed: 212.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam` (239139)
<pre>
</pre>
Command completed. Elapsed time: 0:11:50. Running peak memory: 6.657GB.  
  PID: 239139;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam` (240245)
<pre>
</pre>
Command completed. Elapsed time: 0:11:10. Running peak memory: 6.657GB.  
  PID: 240245;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 11:33:56) elapsed: 1379.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (241320)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 241320;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/plus_TSS.tsv -p ends -c 16 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_plus_TssEnrichment.txt` (241325)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 6.657GB.  
  PID: 241325;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.035GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/minus_TSS.tsv -p ends -c 16 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_minus_TssEnrichment.txt` (241377)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.657GB.  
  PID: 241377;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.25GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_minus_TssEnrichment.txt` (241424)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.657GB.  
  PID: 241424;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_40_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_40_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt` (241449,241450,241451,241452)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 241449;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 241451;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 241450;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 241452;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_keep.txt` (241454)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 241454;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 11:34:39) elapsed: 43.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_ensembl_tss.bed` (241456,241457)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.657GB.  
  PID: 241456;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 241457;	Command: bedtools;	Return code: 0;	Memory used: 0.097GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_ensembl_gene_body.bed` (241461,241462)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 241461;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 241462;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_TSS_density.bed` (241464,241465,241466,241467)
<pre>
</pre>
Command completed. Elapsed time: 0:03:33. Running peak memory: 6.657GB.  
  PID: 241466;	Command: sort;	Return code: 0;	Memory used: 0.012GB  
  PID: 241464;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 241467;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 241465;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_gene_body_density.bed` (241883,241884,241885)
<pre>
</pre>
Command completed. Elapsed time: 0:05:45. Running peak memory: 6.657GB.  
  PID: 241883;	Command: bedtools;	Return code: 0;	Memory used: 0.423GB  
  PID: 241885;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 241884;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/tmppp3luon0` (242501,242502,242503)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 242501;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 242503;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 242502;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/tmppp3luon0 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0999524) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/tmppp3luon0 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_pause_index.bed` (242510)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 242510;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.83	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_pause_index.bed` (242515)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.657GB.  
  PID: 242515;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/K562_PRO-seq_40_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_40_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_pause_index.bed.gz`  

> `pigz -f -p 16 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_pause_index.bed` (242538)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 242538;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 11:44:04) elapsed: 565.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam`
154949505 55262609

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam`
154949505 52574536

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_gene_sort.bed` (243030,243031)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.657GB.  
  PID: 243030;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 243031;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_gene_coverage.bed` (243034)
<pre>
</pre>
Command completed. Elapsed time: 0:05:25. Running peak memory: 6.657GB.  
  PID: 243034;	Command: bedtools;	Return code: 0;	Memory used: 0.421GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/hg38_annotations.bed.gz` (271245)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 271245;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 16 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/hg38_annotations.bed` (271246)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 271246;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 11:54:04) elapsed: 599.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/raw/hg38_annotations.bed` (271256)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.657GB.  
  PID: 271256;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Enhancer_sort.bed` (271258,271259,271260,271261)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.657GB.  
  PID: 271258;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 271259;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 271261;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB  
  PID: 271260;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Enhancer_plus_coverage.bed` (271264)
<pre>
</pre>
Command completed. Elapsed time: 0:01:52. Running peak memory: 6.657GB.  
  PID: 271264;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Enhancer_minus_coverage.bed` (288125)
<pre>
</pre>
Command completed. Elapsed time: 0:01:45. Running peak memory: 6.657GB.  
  PID: 288125;	Command: bedtools;	Return code: 0;	Memory used: 0.032GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_sort.bed` (293753,293754,293755,293756)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 293753;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 293754;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 293756;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 293755;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_plus_coverage.bed` (293758)
<pre>
</pre>
Command completed. Elapsed time: 0:01:47. Running peak memory: 6.657GB.  
  PID: 293758;	Command: bedtools;	Return code: 0;	Memory used: 0.241GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_minus_coverage.bed` (293904)
<pre>
</pre>
Command completed. Elapsed time: 0:01:47. Running peak memory: 6.657GB.  
  PID: 293904;	Command: bedtools;	Return code: 0;	Memory used: 0.193GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_Flanking_Region"` (313533)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 313533;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed` (313564,313565,313566,313568)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.657GB.  
  PID: 313564;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 313566;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 313565;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 313568;	Command: bedtools;	Return code: 0;	Memory used: 0.051GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_Flanking_Region_plus_coverage.bed` (313901)
<pre>
</pre>
Command completed. Elapsed time: 0:01:45. Running peak memory: 6.657GB.  
  PID: 313901;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_Flanking_Region_minus_coverage.bed` (338297)
<pre>
</pre>
Command completed. Elapsed time: 0:01:39. Running peak memory: 6.657GB.  
  PID: 338297;	Command: bedtools;	Return code: 0;	Memory used: 0.069GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/5_UTR"` (338439)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 338439;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/5_UTR_sort.bed` (338440,338441,338442,338443)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.657GB.  
  PID: 338440;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 338441;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 338443;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 338442;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_5_UTR_plus_coverage.bed` (338447)
<pre>
</pre>
Command completed. Elapsed time: 0:01:40. Running peak memory: 6.657GB.  
  PID: 338447;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_5_UTR_minus_coverage.bed` (338725)
<pre>
</pre>
Command completed. Elapsed time: 0:01:40. Running peak memory: 6.657GB.  
  PID: 338725;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/3_UTR"` (369449)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 369449;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/3_UTR_sort.bed` (369517,369530,369533,369546)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.657GB.  
  PID: 369517;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 369530;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 369546;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 369533;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_3_UTR_plus_coverage.bed` (370164)
<pre>
</pre>
Command completed. Elapsed time: 0:01:45. Running peak memory: 6.657GB.  
  PID: 370164;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_3_UTR_minus_coverage.bed` (382951)
<pre>
</pre>
Command completed. Elapsed time: 0:01:46. Running peak memory: 6.657GB.  
  PID: 382951;	Command: bedtools;	Return code: 0;	Memory used: 0.069GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Exon_sort.bed` (383374,383375,383376,383377)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.657GB.  
  PID: 383374;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 383375;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 383377;	Command: bedtools;	Return code: 0;	Memory used: 0.172GB  
  PID: 383376;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Exon_plus_coverage.bed` (383381)
<pre>
</pre>
Command completed. Elapsed time: 0:01:54. Running peak memory: 6.657GB.  
  PID: 383381;	Command: bedtools;	Return code: 0;	Memory used: 0.194GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Exon_minus_coverage.bed` (383482)
<pre>
</pre>
Command completed. Elapsed time: 0:01:53. Running peak memory: 6.657GB.  
  PID: 383482;	Command: bedtools;	Return code: 0;	Memory used: 0.081GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Intron_sort.bed` (427749,427750,427751,427752)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.657GB.  
  PID: 427749;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 427751;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 427750;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 427752;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Intron_plus_coverage.bed` (427756)
<pre>
</pre>
Command completed. Elapsed time: 0:02:06. Running peak memory: 6.657GB.  
  PID: 427756;	Command: bedtools;	Return code: 0;	Memory used: 0.188GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Intron_minus_coverage.bed` (427932)
<pre>
</pre>
Command completed. Elapsed time: 0:02:04. Running peak memory: 6.657GB.  
  PID: 427932;	Command: bedtools;	Return code: 0;	Memory used: 0.181GB


### Plot cFRiF/FRiF (06-15 12:19:37) elapsed: 1534.0 _TIME_


> `samtools view -@ 16 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_40 -z 3099922541 -n 77100431 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Intron_plus_coverage.bed` (447316)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 6.657GB.  
  PID: 447316;	Command: Rscript;	Return code: 0;	Memory used: 0.436GB

> `cFRiF`	QC_hg38/K562_PRO-seq_40_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_40_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_40 -z 3099922541 -n 77100431 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_Intron_plus_coverage.bed` (2927)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.657GB.  
  PID: 2927;	Command: Rscript;	Return code: 0;	Memory used: 0.462GB

> `FRiF`	QC_hg38/K562_PRO-seq_40_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_40_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 12:20:56) elapsed: 78.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_exons_sort.bed` (13858,13862)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.657GB.  
  PID: 13862;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB  
  PID: 13858;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_introns_sort.bed` (14036,14037,14038)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.657GB.  
  PID: 14036;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 14038;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 14037;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exons_coverage.bed` (14045)
<pre>
</pre>
Command completed. Elapsed time: 0:03:38. Running peak memory: 6.657GB.  
  PID: 14045;	Command: bedtools;	Return code: 0;	Memory used: 0.121GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_introns_coverage.bed` (35913)
<pre>
</pre>
Command completed. Elapsed time: 0:04:17. Running peak memory: 6.657GB.  
  PID: 35913;	Command: bedtools;	Return code: 0;	Memory used: 0.198GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/154.949505)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exons_rpkm.bed` (58844,58845,58846)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.657GB.  
  PID: 58844;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 58846;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 58845;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/154.949505)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_introns_rpkm.bed` (58848,58849,58850)
<pre>
psutil.NoSuchProcess process no longer exists (pid=58851)
Warning: couldn't add memory use for process: 58850
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.657GB.  
  PID: 58848;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 58850;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 58849;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exon_intron_ratios.bed` (58852,58853,58854)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 58852;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 58854;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 58853;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.38	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exon_intron_ratios.bed --annotate` (58861)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.657GB.  
  PID: 58861;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_40_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_40_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exon_intron_ratios.bed.gz`  

> `pigz -f -p 16 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/QC_hg38/K562_PRO-seq_40_exon_intron_ratios.bed` (58909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.657GB.  
  PID: 58909;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (06-15 12:29:07) elapsed: 492.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam` (58918)
<pre>
</pre>
Command completed. Elapsed time: 0:01:13. Running peak memory: 6.657GB.  
  PID: 58918;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_plus_smooth_body_0-mer.bw -p 10 --variable-step --tail-edge --scale 154949505.0` (93416)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_40_plus_cuttrace_2if4cqjs'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 94 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 101 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 101 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_plus_exact_body_0-mer.bw'
Merging 101 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:07. Running peak memory: 6.657GB.  
  PID: 93416;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.463GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam` (193278)
<pre>
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 6.657GB.  
  PID: 193278;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_minus_smooth_body_0-mer.bw -p 10 --variable-step --tail-edge --scale 154949505.0` (193351)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/aligned_hg38/K562_PRO-seq_40_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_40_minus_cuttrace_ze8wk8bh'
Processing with 5 cores...
stdin is empty of data
Discarding 102 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 93 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 93 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_minus_exact_body_0-mer.bw'
Merging 93 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_40/signal_hg38/K562_PRO-seq_40_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:22. Running peak memory: 6.657GB.  
  PID: 193351;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.498GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  5:38:58
*  Total elapsed time (all runs):  7:10:10
*         Peak memory (this run):  6.6573 GB
*        Pipeline completed time: 2020-06-15 12:49:59
