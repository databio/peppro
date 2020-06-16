### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_30 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_30pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 16 -M 16000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-12c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/
*  Pipeline started at:   (06-15 07:11:00) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_30pct.fastq.gz']`
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
*        `sample_name`:  `K562_PRO-seq_30`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_30pct.fastq.gz

> `File_mb`	10969.75	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:11:01) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/K562_PRO-seq_30.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_PRO_30pct.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/K562_PRO-seq_30.fastq.gz` (297106)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 297106;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/K562_PRO-seq_30.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1.fastq`  

> `pigz -f -p 16 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/K562_PRO-seq_30.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1.fastq` (297107)
<pre>
</pre>
Command completed. Elapsed time: 0:03:54. Running peak memory: 0.003GB.  
  PID: 297107;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	148988721	PEPPRO	_RES_

> `Fastq_reads`	148988721	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/K562_PRO-seq_30.fastq.gz']

### FASTQ processing:  (06-15 07:18:41) elapsed: 461.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1_processed.fastq`  

> `(cutadapt -j 16 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt` (297777)
<pre>
</pre>
Command completed. Elapsed time: 0:05:08. Running peak memory: 6.459GB.  
  PID: 297777;	Command: cutadapt;	Return code: 0;	Memory used: 6.459GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1_processed.fastq` (298327,298328)
<pre>
</pre>
Command completed. Elapsed time: 0:03:54. Running peak memory: 6.459GB.  
  PID: 298327;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 298328;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	126928151.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	3729909.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.5035	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastqc/K562_PRO-seq_30_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (299217)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 299217;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	145258812	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1_processed.fastq` (299268)
<pre>
Started analysis of K562_PRO-seq_30_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_30_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_30_R1_processed.fastq
Analysis complete for K562_PRO-seq_30_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:07:37. Running peak memory: 6.459GB.  
  PID: 299268;	Command: fastqc;	Return code: 0;	Memory used: 0.227GB

> `FastQC report r1`	fastqc/K562_PRO-seq_30_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/processed_R1.flag` (299933)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 299933;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-15 07:39:17) elapsed: 1236.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt` (299940)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.459GB.  
  PID: 299940;	Command: Rscript;	Return code: 0;	Memory used: 0.112GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_30_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_30_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:39:25) elapsed: 8.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/cutadapt/K562_PRO-seq_30_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2312	PEPPRO	_RES_

### Prealignments (06-15 07:39:25) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:39:25) elapsed: 0.0 _TIME_


> `(bowtie2 -p 16 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq_30 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/fastq/K562_PRO-seq_30_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/prealignments/K562_PRO-seq_30_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
145258812 reads; of these:
  145258812 (100.00%) were unpaired; of these:
    131853845 (90.77%) aligned 0 times
    13404967 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	13404967.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-15 07:58:31) elapsed: 1145.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam`  

> `bowtie2 -p 16 --very-sensitive --rg-id K562_PRO-seq_30 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/prealignments/K562_PRO-seq_30_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/tmpwnwe2w9n -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_temp.bam` (302092,302098,302099)
<pre>
131853845 reads; of these:
  131853845 (100.00%) were unpaired; of these:
    1619546 (1.23%) aligned 0 times
    98460662 (74.67%) aligned exactly 1 time
    31773637 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 42 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 1:17:27. Running peak memory: 6.459GB.  
  PID: 302098;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 302092;	Command: bowtie2;	Return code: 0;	Memory used: 3.817GB  
  PID: 302099;	Command: samtools;	Return code: 0;	Memory used: 0.96GB


> `samtools view -q 10 -b -@ 16 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam` (311201)
<pre>
</pre>
Command completed. Elapsed time: 0:03:49. Running peak memory: 6.459GB.  
  PID: 311201;	Command: samtools;	Return code: 0;	Memory used: 0.022GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	130234299	PEPPRO	_RES_

> `QC_filtered_reads`	14019843	PEPPRO	_RES_

> `Aligned_reads`	116214456	PEPPRO	_RES_

> `Alignment_rate`	80.01	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	10.45	PEPPRO	_RES_

### Compress all unmapped read files (06-15 09:52:13) elapsed: 6823.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/prealignments/K562_PRO-seq_30_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/prealignments/K562_PRO-seq_30_human_rDNA_unmap.fq` (315212)
<pre>
</pre>
Command completed. Elapsed time: 0:03:26. Running peak memory: 6.459GB.  
  PID: 315212;	Command: pigz;	Return code: 0;	Memory used: 0.012GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_temp.bam` (315617)
<pre>
</pre>
Command completed. Elapsed time: 0:02:13. Running peak memory: 6.459GB.  
  PID: 315617;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	2744825	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam` (315732)
<pre>
</pre>
Command completed. Elapsed time: 0:01:56. Running peak memory: 6.459GB.  
  PID: 315732;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/chr_sizes.bed` (315839,315840,315841,315842)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 315841;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 315839;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 315842;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 315840;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/chr_sizes.bed -b -@ 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_noMT.bam` (315844)
<pre>
</pre>
Command completed. Elapsed time: 0:02:26. Running peak memory: 6.459GB.  
  PID: 315844;	Command: samtools;	Return code: 0;	Memory used: 0.023GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam` (316322)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 316322;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam` (316324)
<pre>
</pre>
Command completed. Elapsed time: 0:01:54. Running peak memory: 6.459GB.  
  PID: 316324;	Command: samtools;	Return code: 0;	Memory used: 0.019GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 10:10:33) elapsed: 1100.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam -c 16 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_bamQC.tsv` (317358)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/tmp_K562_PRO-seq_30_sort_qxylztwb'
Processing with 16 cores...
Discarding 87 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 108 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:02:12. Running peak memory: 6.459GB.  
  PID: 317358;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 2.928GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_bamQC.tsv`

> `NRF`	0.6	PEPPRO	_RES_

> `PBC1`	0.79	PEPPRO	_RES_

> `PBC2`	6.54	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_unmap.bam`  

> `samtools view -b -@ 16 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_unmap.bam` (317628)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.459GB.  
  PID: 317628;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools view -c -f 4 -@ 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_temp.bam`

> `Unmapped_reads`	1619546	PEPPRO	_RES_

### Split BAM by strand (06-15 10:13:33) elapsed: 179.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam` (317704)
<pre>
</pre>
Command completed. Elapsed time: 0:09:25. Running peak memory: 6.459GB.  
  PID: 317704;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam` (318739)
<pre>
</pre>
Command completed. Elapsed time: 0:09:12. Running peak memory: 6.459GB.  
  PID: 318739;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 10:32:10) elapsed: 1118.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (319770)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 319770;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/plus_TSS.tsv -p ends -c 16 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_plus_TssEnrichment.txt` (319775)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.459GB.  
  PID: 319775;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.237GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/minus_TSS.tsv -p ends -c 16 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_minus_TssEnrichment.txt` (319829)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.459GB.  
  PID: 319829;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.737GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_minus_TssEnrichment.txt` (319878)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.459GB.  
  PID: 319878;	Command: Rscript;	Return code: 0;	Memory used: 0.285GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_30_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_30_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt` (319904,319905,319906,319907)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 319904;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 319906;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 319905;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 319907;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_keep.txt` (319909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 319909;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 10:32:46) elapsed: 36.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_ensembl_tss.bed` (319911,319912)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.459GB.  
  PID: 319911;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 319912;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_ensembl_gene_body.bed` (319917,319918)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 319917;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 319918;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_TSS_density.bed` (319920,319921,319922,319923)
<pre>
</pre>
Command completed. Elapsed time: 0:03:00. Running peak memory: 6.459GB.  
  PID: 319920;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 319922;	Command: sort;	Return code: 0;	Memory used: 0.01GB  
  PID: 319921;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 319923;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_gene_body_density.bed` (320340,320341,320342)
<pre>
</pre>
Command completed. Elapsed time: 0:05:23. Running peak memory: 6.459GB.  
  PID: 320342;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 320340;	Command: bedtools;	Return code: 0;	Memory used: 0.322GB  
  PID: 320341;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/tmpm9nnjnm6` (320924,320925,320926)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 320924;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 320926;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 320925;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/tmpm9nnjnm6 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0790864) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/tmpm9nnjnm6 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_pause_index.bed` (320932)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 320932;	Command: awk;	Return code: 0;	Memory used: 0.004GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.73	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_pause_index.bed` (320937)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.459GB.  
  PID: 320937;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/K562_PRO-seq_30_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_30_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_pause_index.bed.gz`  

> `pigz -f -p 16 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_pause_index.bed` (320961)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 320961;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 10:41:17) elapsed: 511.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam`
116214456 41446005

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam`
116214456 39433548

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_gene_sort.bed` (321180,321181)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 321180;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 321181;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_gene_coverage.bed` (321184)
<pre>
</pre>
Command completed. Elapsed time: 0:05:09. Running peak memory: 6.459GB.  
  PID: 321184;	Command: bedtools;	Return code: 0;	Memory used: 0.321GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/hg38_annotations.bed.gz` (321712)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 321712;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 16 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/hg38_annotations.bed` (321713)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 321713;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 10:49:57) elapsed: 520.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/raw/hg38_annotations.bed` (321722)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 321722;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Enhancer_sort.bed` (321724,321725,321726,321727)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 321724;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 321725;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 321727;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 321726;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Enhancer_plus_coverage.bed` (321729)
<pre>
</pre>
Command completed. Elapsed time: 0:01:25. Running peak memory: 6.459GB.  
  PID: 321729;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Enhancer_minus_coverage.bed` (322096)
<pre>
</pre>
Command completed. Elapsed time: 0:01:30. Running peak memory: 6.459GB.  
  PID: 322096;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_sort.bed` (322175,322176,322177,322178)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 322175;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 322176;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 322178;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 322177;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_plus_coverage.bed` (322180)
<pre>
</pre>
Command completed. Elapsed time: 0:01:48. Running peak memory: 6.459GB.  
  PID: 322180;	Command: bedtools;	Return code: 0;	Memory used: 0.182GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_minus_coverage.bed` (322290)
<pre>
</pre>
Command completed. Elapsed time: 0:01:40. Running peak memory: 6.459GB.  
  PID: 322290;	Command: bedtools;	Return code: 0;	Memory used: 0.145GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_Flanking_Region"` (322613)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 322613;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed` (322614,322615,322616,322618)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 322614;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 322616;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 322615;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 322618;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_Flanking_Region_plus_coverage.bed` (322622)
<pre>
</pre>
Command completed. Elapsed time: 0:01:29. Running peak memory: 6.459GB.  
  PID: 322622;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_Flanking_Region_minus_coverage.bed` (322804)
<pre>
</pre>
Command completed. Elapsed time: 0:01:39. Running peak memory: 6.459GB.  
  PID: 322804;	Command: bedtools;	Return code: 0;	Memory used: 0.054GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/5_UTR"` (322918)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 322918;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/5_UTR_sort.bed` (322919,322920,322921,322922)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 322919;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 322920;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 322922;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 322921;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_5_UTR_plus_coverage.bed` (322925)
<pre>
</pre>
Command completed. Elapsed time: 0:01:32. Running peak memory: 6.459GB.  
  PID: 322925;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_5_UTR_minus_coverage.bed` (323331)
<pre>
</pre>
Command completed. Elapsed time: 0:01:28. Running peak memory: 6.459GB.  
  PID: 323331;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/3_UTR"` (323429)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 323429;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/3_UTR_sort.bed` (323430,323431,323432,323434)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 323430;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 323431;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 323434;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 323432;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_3_UTR_plus_coverage.bed` (323437)
<pre>
</pre>
Command completed. Elapsed time: 0:01:26. Running peak memory: 6.459GB.  
  PID: 323437;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_3_UTR_minus_coverage.bed` (323563)
<pre>
</pre>
Command completed. Elapsed time: 0:01:24. Running peak memory: 6.459GB.  
  PID: 323563;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Exon_sort.bed` (323851,323852,323853,323854)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.459GB.  
  PID: 323851;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 323852;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 323854;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 323853;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Exon_plus_coverage.bed` (323858)
<pre>
</pre>
Command completed. Elapsed time: 0:01:44. Running peak memory: 6.459GB.  
  PID: 323858;	Command: bedtools;	Return code: 0;	Memory used: 0.149GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Exon_minus_coverage.bed` (323949)
<pre>
</pre>
Command completed. Elapsed time: 0:01:43. Running peak memory: 6.459GB.  
  PID: 323949;	Command: bedtools;	Return code: 0;	Memory used: 0.043GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Intron_sort.bed` (324042,324043,324044,324045)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.459GB.  
  PID: 324042;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 324044;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 324043;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 324045;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Intron_plus_coverage.bed` (324049)
<pre>
</pre>
Command completed. Elapsed time: 0:02:10. Running peak memory: 6.459GB.  
  PID: 324049;	Command: bedtools;	Return code: 0;	Memory used: 0.145GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Intron_minus_coverage.bed` (324452)
<pre>
</pre>
Command completed. Elapsed time: 0:02:09. Running peak memory: 6.459GB.  
  PID: 324452;	Command: bedtools;	Return code: 0;	Memory used: 0.139GB


### Plot cFRiF/FRiF (06-15 11:13:15) elapsed: 1398.0 _TIME_


> `samtools view -@ 16 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_30 -z 3099922541 -n 57825179 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Intron_plus_coverage.bed` (324595)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 6.459GB.  
  PID: 324595;	Command: Rscript;	Return code: 0;	Memory used: 0.849GB

> `cFRiF`	QC_hg38/K562_PRO-seq_30_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_30_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_30 -z 3099922541 -n 57825179 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_Intron_plus_coverage.bed` (324652)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.459GB.  
  PID: 324652;	Command: Rscript;	Return code: 0;	Memory used: 0.432GB

> `FRiF`	QC_hg38/K562_PRO-seq_30_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_30_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 11:14:43) elapsed: 88.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_exons_sort.bed` (324699,324700)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.459GB.  
  PID: 324700;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 324699;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_introns_sort.bed` (324708,324709,324710)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.459GB.  
  PID: 324708;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 324710;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 324709;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exons_coverage.bed` (324717)
<pre>
</pre>
Command completed. Elapsed time: 0:03:01. Running peak memory: 6.459GB.  
  PID: 324717;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_introns_coverage.bed` (325110)
<pre>
</pre>
Command completed. Elapsed time: 0:04:23. Running peak memory: 6.459GB.  
  PID: 325110;	Command: bedtools;	Return code: 0;	Memory used: 0.152GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/116.214456)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exons_rpkm.bed` (325651,325652,325653)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 325651;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 325653;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 325652;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/116.214456)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_introns_rpkm.bed` (325656,325657,325658)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.459GB.  
  PID: 325656;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 325658;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 325657;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exon_intron_ratios.bed` (325663,325664,325665)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 325663;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 325665;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 325664;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.38	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exon_intron_ratios.bed --annotate` (325671)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.459GB.  
  PID: 325671;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_30_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_30_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exon_intron_ratios.bed.gz`  

> `pigz -f -p 16 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/QC_hg38/K562_PRO-seq_30_exon_intron_ratios.bed` (325694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.459GB.  
  PID: 325694;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 11:22:24) elapsed: 461.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam` (325704)
<pre>
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 6.459GB.  
  PID: 325704;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_plus_smooth_body_0-mer.bw -p 10 --variable-step --tail-edge --scale 116214456.0` (325757)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_30_plus_cuttrace_y6jjp0nj'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 94 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 101 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 101 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_plus_exact_body_0-mer.bw'
Merging 101 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:22. Running peak memory: 6.459GB.  
  PID: 325757;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.804GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam` (327849)
<pre>
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 6.459GB.  
  PID: 327849;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_minus_smooth_body_0-mer.bw -p 10 --variable-step --tail-edge --scale 116214456.0` (327914)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/aligned_hg38/K562_PRO-seq_30_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_30_minus_cuttrace_n5pv2oxe'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
Discarding 102 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 93 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 93 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_minus_exact_body_0-mer.bw'
Merging 93 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_30/signal_hg38/K562_PRO-seq_30_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:34. Running peak memory: 6.459GB.  
  PID: 327914;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.761GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  4:32:19
*  Total elapsed time (all runs):  5:43:26
*         Peak memory (this run):  6.459 GB
*        Pipeline completed time: 2020-06-15 11:43:19
