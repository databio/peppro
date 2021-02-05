### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_70 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_70pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-23a
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/
*  Pipeline started at:   (06-15 07:17:11) elapsed: 0.0 _TIME_

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
*              `cores`:  `12`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_70pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_70`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_70pctRNA.fastq.gz

> `File_mb`	5404.02	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:12) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/K562_RNA-seq_70.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_70pctRNA.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/K562_RNA-seq_70.fastq.gz` (13570)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 13570;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/K562_RNA-seq_70.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/K562_RNA-seq_70.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1.fastq` (13571)
<pre>
</pre>
Command completed. Elapsed time: 0:04:42. Running peak memory: 0.002GB.  
  PID: 13571;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	70000000	PEPPRO	_RES_

> `Fastq_reads`	70000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/K562_RNA-seq_70.fastq.gz']

### FASTQ processing:  (06-15 07:23:53) elapsed: 401.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt` (14198)
<pre>
</pre>
Command completed. Elapsed time: 0:02:56. Running peak memory: 3.999GB.  
  PID: 14198;	Command: cutadapt;	Return code: 0;	Memory used: 3.999GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_processed.fastq` (14579,14580)
<pre>
</pre>
Command completed. Elapsed time: 0:02:06. Running peak memory: 3.999GB.  
  PID: 14579;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 14580;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	35291630.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	524857.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	0.7498	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastqc/K562_RNA-seq_70_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (15064)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 15064;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	69475143	PEPPRO	_RES_

> `Trim_loss_rate`	0.75	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_processed.fastq` (15098)
<pre>
Started analysis of K562_RNA-seq_70_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_70_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_70_R1_processed.fastq
Analysis complete for K562_RNA-seq_70_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:04:19. Running peak memory: 3.999GB.  
  PID: 15098;	Command: fastqc;	Return code: 0;	Memory used: 0.244GB

> `FastQC report r1`	fastqc/K562_RNA-seq_70_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/processed_R1.flag` (15668)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 15668;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Plot adapter insertion distribution (06-15 07:35:30) elapsed: 697.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt` (15669)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.999GB.  
  PID: 15669;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_70_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_70_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:35:35) elapsed: 5.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/cutadapt/K562_RNA-seq_70_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2308	PEPPRO	_RES_

### Prealignments (06-15 07:35:35) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:35:35) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_70 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/fastq/K562_RNA-seq_70_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/prealignments/K562_RNA-seq_70_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
69475143 reads; of these:
  69475143 (100.00%) were unpaired; of these:
    66449603 (95.65%) aligned 0 times
    3025540 (4.35%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
4.35% overall alignment rate

> `Aligned_reads_human_rDNA`	3025540.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	4.35	PEPPRO	_RES_

### Map to genome (06-15 07:44:33) elapsed: 538.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_RNA-seq_70 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/prealignments/K562_RNA-seq_70_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/tmpzz7s0tum -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam` (16477,16483,16484)
<pre>
66449603 reads; of these:
  66449603 (100.00%) were unpaired; of these:
    8296734 (12.49%) aligned 0 times
    37395744 (56.28%) aligned exactly 1 time
    20757125 (31.24%) aligned >1 times
87.51% overall alignment rate
[bam_sort_core] merging from 21 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:35:47. Running peak memory: 3.999GB.  
  PID: 16483;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 16477;	Command: bowtie2;	Return code: 0;	Memory used: 3.713GB  
  PID: 16484;	Command: samtools;	Return code: 0;	Memory used: 0.889GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam` (20398)
<pre>
</pre>
Command completed. Elapsed time: 0:03:15. Running peak memory: 3.999GB.  
  PID: 20398;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	58152869	PEPPRO	_RES_

> `QC_filtered_reads`	10854484	PEPPRO	_RES_

> `Aligned_reads`	47298385	PEPPRO	_RES_

> `Alignment_rate`	68.08	PEPPRO	_RES_

> `Total_efficiency`	67.57	PEPPRO	_RES_

> `Read_depth`	7.29	PEPPRO	_RES_

### Compress all unmapped read files (06-15 08:41:55) elapsed: 3441.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/prealignments/K562_RNA-seq_70_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/prealignments/K562_RNA-seq_70_human_rDNA_unmap.fq` (22407)
<pre>
</pre>
Command completed. Elapsed time: 0:02:45. Running peak memory: 3.999GB.  
  PID: 22407;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam` (22566)
<pre>
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 3.999GB.  
  PID: 22566;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	3321940	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam` (22823)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 3.999GB.  
  PID: 22823;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/chr_sizes.bed` (22868,22869,22870,22871)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 22868;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 22870;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 22869;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 22871;	Command: grep;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_noMT.bam` (22873)
<pre>
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 3.999GB.  
  PID: 22873;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam` (22934)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 22934;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam` (22935)
<pre>
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 3.999GB.  
  PID: 22935;	Command: samtools;	Return code: 0;	Memory used: 0.016GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 08:50:51) elapsed: 537.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_bamQC.tsv` (23334)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/tmp_K562_RNA-seq_70_sort_lm_6ny4f'
Processing with 12 cores...
Discarding 97 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 98 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:01:04. Running peak memory: 3.999GB.  
  PID: 23334;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.396GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_bamQC.tsv`

> `NRF`	0.65	PEPPRO	_RES_

> `PBC1`	0.88	PEPPRO	_RES_

> `PBC2`	13.92	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_unmap.bam` (23424)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.999GB.  
  PID: 23424;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_temp.bam`

> `Unmapped_reads`	8296734	PEPPRO	_RES_

### Split BAM by strand (06-15 08:52:26) elapsed: 94.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam` (23503)
<pre>
</pre>
Command completed. Elapsed time: 0:03:42. Running peak memory: 3.999GB.  
  PID: 23503;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam` (23880)
<pre>
</pre>
Command completed. Elapsed time: 0:03:37. Running peak memory: 3.999GB.  
  PID: 23880;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 08:59:45) elapsed: 439.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (24065)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 24065;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_plus_TssEnrichment.txt` (24066)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.999GB.  
  PID: 24066;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.519GB


> `TSS_coding_score`	17.7	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_minus_TssEnrichment.txt` (24101)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.999GB.  
  PID: 24101;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.581GB


> `TSS_non-coding_score`	3.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_minus_TssEnrichment.txt` (24172)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.999GB.  
  PID: 24172;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_70_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_70_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt` (24404,24405,24406,24407)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 24404;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 24406;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 24405;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 24407;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt` (24410)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 24410;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 09:00:09) elapsed: 24.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_tss.bed` (24412,24413)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.999GB.  
  PID: 24412;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 24413;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_gene_body.bed` (24417,24418)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 24417;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 24418;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_TSS_density.bed` (24420,24421,24422,24423)
<pre>
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 3.999GB.  
  PID: 24421;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 24423;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 24420;	Command: bedtools;	Return code: 0;	Memory used: 0.053GB  
  PID: 24422;	Command: sort;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_gene_body_density.bed` (24497,24498,24499)
<pre>
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 3.999GB.  
  PID: 24497;	Command: bedtools;	Return code: 0;	Memory used: 0.244GB  
  PID: 24499;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 24498;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/tmp43ogt6mm` (24588,24589,24590)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 24588;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 24590;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 24589;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/tmp43ogt6mm | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0384066) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/tmp43ogt6mm > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed` (24596)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 24596;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	7.99	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed` (24601)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.999GB.  
  PID: 24601;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/K562_RNA-seq_70_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_70_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_pause_index.bed` (24622)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 24622;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 09:03:10) elapsed: 181.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam`
47298385 18212495

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam`
47298385 18149986

> `Minus_FRiP`	0.38	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_gene_sort.bed` (24740,24741)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.999GB.  
  PID: 24740;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 24741;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_gene_coverage.bed` (24743)
<pre>
</pre>
Command completed. Elapsed time: 0:01:39. Running peak memory: 3.999GB.  
  PID: 24743;	Command: bedtools;	Return code: 0;	Memory used: 0.198GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed.gz` (25016)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 25016;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed` (25017)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 25017;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 09:06:09) elapsed: 179.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/raw/hg38_annotations.bed` (25026)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.999GB.  
  PID: 25026;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer_sort.bed` (25028,25029,25030,25031)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.999GB.  
  PID: 25028;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 25029;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 25031;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 25030;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_plus_coverage.bed` (25034)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.999GB.  
  PID: 25034;	Command: bedtools;	Return code: 0;	Memory used: 0.028GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_minus_coverage.bed` (25064)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.999GB.  
  PID: 25064;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_sort.bed` (25096,25097,25098,25099)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 25096;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 25097;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 25099;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 25098;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_plus_coverage.bed` (25101)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.999GB.  
  PID: 25101;	Command: bedtools;	Return code: 0;	Memory used: 0.141GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_minus_coverage.bed` (25134)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 3.999GB.  
  PID: 25134;	Command: bedtools;	Return code: 0;	Memory used: 0.07GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region"` (25169)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 25169;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed` (25170,25171,25172,25173)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.999GB.  
  PID: 25170;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 25172;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 25171;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 25173;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_plus_coverage.bed` (25176)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 3.999GB.  
  PID: 25176;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_minus_coverage.bed` (25209)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.999GB.  
  PID: 25209;	Command: bedtools;	Return code: 0;	Memory used: 0.051GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR"` (25241)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 25241;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR_sort.bed` (25242,25243,25244,25245)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.999GB.  
  PID: 25242;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 25243;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 25245;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 25244;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_plus_coverage.bed` (25247)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.999GB.  
  PID: 25247;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_minus_coverage.bed` (25515)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.999GB.  
  PID: 25515;	Command: bedtools;	Return code: 0;	Memory used: 0.029GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR"` (25547)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 25547;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR_sort.bed` (25548,25549,25550,25551)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.999GB.  
  PID: 25548;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 25549;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 25551;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 25550;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_plus_coverage.bed` (25554)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 3.999GB.  
  PID: 25554;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_minus_coverage.bed` (25587)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 3.999GB.  
  PID: 25587;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon_sort.bed` (25620,25622,25623,25624)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.999GB.  
  PID: 25620;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 25622;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 25624;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 25623;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_plus_coverage.bed` (25628)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 3.999GB.  
  PID: 25628;	Command: bedtools;	Return code: 0;	Memory used: 0.177GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_minus_coverage.bed` (25669)
<pre>
</pre>
Command completed. Elapsed time: 0:00:43. Running peak memory: 3.999GB.  
  PID: 25669;	Command: bedtools;	Return code: 0;	Memory used: 0.063GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron_sort.bed` (25707,25708,25709,25710)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.999GB.  
  PID: 25707;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 25709;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 25708;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 25710;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_plus_coverage.bed` (25713)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 3.999GB.  
  PID: 25713;	Command: bedtools;	Return code: 0;	Memory used: 0.132GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_minus_coverage.bed` (25752)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 3.999GB.  
  PID: 25752;	Command: bedtools;	Return code: 0;	Memory used: 0.067GB


### Plot cFRiF/FRiF (06-15 09:15:05) elapsed: 537.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_70 -z 3099922541 -n 22361014 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_plus_coverage.bed` (25989)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.999GB.  
  PID: 25989;	Command: Rscript;	Return code: 0;	Memory used: 0.445GB

> `cFRiF`	QC_hg38/K562_RNA-seq_70_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_70_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_70 -z 3099922541 -n 22361014 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_Intron_plus_coverage.bed` (26039)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 3.999GB.  
  PID: 26039;	Command: Rscript;	Return code: 0;	Memory used: 0.445GB

> `FRiF`	QC_hg38/K562_RNA-seq_70_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_70_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 09:16:15) elapsed: 70.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_exons_sort.bed` (26073,26074)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.999GB.  
  PID: 26074;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB  
  PID: 26073;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_introns_sort.bed` (26080,26081,26082)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.999GB.  
  PID: 26080;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 26082;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 26081;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_coverage.bed` (26088)
<pre>
</pre>
Command completed. Elapsed time: 0:01:22. Running peak memory: 3.999GB.  
  PID: 26088;	Command: bedtools;	Return code: 0;	Memory used: 0.162GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_coverage.bed` (26160)
<pre>
</pre>
Command completed. Elapsed time: 0:01:22. Running peak memory: 3.999GB.  
  PID: 26160;	Command: bedtools;	Return code: 0;	Memory used: 0.151GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/47.298385)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_rpkm.bed` (26232,26233,26234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.999GB.  
  PID: 26232;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 26234;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 26233;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/47.298385)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_rpkm.bed` (26236,26237,26238)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.999GB.  
  PID: 26236;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 26238;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 26237;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed` (26241,26242,26243)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 26241;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 26243;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 26242;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	5.92	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed --annotate` (26249)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.999GB.  
  PID: 26249;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_70_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_70_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/QC_hg38/K562_RNA-seq_70_exon_intron_ratios.bed` (26271)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.999GB.  
  PID: 26271;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 09:19:16) elapsed: 181.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam` (26279)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.999GB.  
  PID: 26279;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 47298385.0` (26299)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_70_plus_cuttrace_vp0t6xfp'
Processing with 4 cores...
Discarding 108 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 87 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 87 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_exact_body_0-mer.bw'
Merging 87 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:00. Running peak memory: 3.999GB.  
  PID: 26299;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.59GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam` (27989)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.999GB.  
  PID: 27989;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 47298385.0` (28010)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/aligned_hg38/K562_RNA-seq_70_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_70_minus_cuttrace_8i5_a9as'
Processing with 4 cores...
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_minus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_70/signal_hg38/K562_RNA-seq_70_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:57. Running peak memory: 3.999GB.  
  PID: 28010;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.608GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:18:47
*  Total elapsed time (all runs):  2:45:24
*         Peak memory (this run):  3.9988 GB
*        Pipeline completed time: 2020-06-15 09:35:58
