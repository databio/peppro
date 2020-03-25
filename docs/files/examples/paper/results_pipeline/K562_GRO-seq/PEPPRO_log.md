### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_GRO-seq --genome hg38 --input /project/shefflab/data/sra_fastq/SRR1552484.fastq.gz --single-or-paired single --protocol GRO --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-ba26-16
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/
*  Pipeline started at:   (02-27 09:24:34) elapsed: 0.0 _TIME_

### Version log:

*       Python version:  3.6.5
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.1
*        Pipeline hash:  2fe0657f50e41000560af043f4914b3a240296f2
*      Pipeline branch:  * dev
*        Pipeline date:  2020-02-27 09:20:39 -0500

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `cutadapt`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `8`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data/sra_fastq/SRR1552484.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `16000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed/peppro/paper/dev4/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `GRO`
*            `recover`:  `False`
*        `sample_name`:  `K562_GRO-seq`
*              `scale`:  `False`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `single`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `0`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /project/shefflab/data/sra_fastq/SRR1552484.fastq.gz

> `File_mb`	1322.25	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected GRO input

### Merge/link and fastq conversion:  (02-27 09:24:35) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz`  

> `ln -sf /project/shefflab/data/sra_fastq/SRR1552484.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz` (96232)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 96232;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1.fastq` (96233)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 0.003GB.  
  PID: 96233;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	30363736	PEPPRO	_RES_

> `Fastq_reads`	30363736	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/K562_GRO-seq.fastq.gz']

### FASTQ processing:  (02-27 09:25:40) elapsed: 66.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_processed.fastq`  

> `(cutadapt -j 8 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt` (96510)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 0.28GB.  
  PID: 96510;	Command: cutadapt;	Return code: 0;	Memory used: 0.28GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_noadap.fastq | seqtk seq -L 2 - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_processed.fastq` (96573,96574)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 0.28GB.  
  PID: 96573;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 96574;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	22296806.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	1041369.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	3.4296	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastqc/K562_GRO-seq_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (96626)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.28GB.  
  PID: 96626;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	29322367	PEPPRO	_RES_

> `Trim_loss_rate`	3.43	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_processed.fastq` (96637)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of K562_GRO-seq_R1_processed.fastq
Approx 5% complete for K562_GRO-seq_R1_processed.fastq
Approx 10% complete for K562_GRO-seq_R1_processed.fastq
Approx 15% complete for K562_GRO-seq_R1_processed.fastq
Approx 20% complete for K562_GRO-seq_R1_processed.fastq
Approx 25% complete for K562_GRO-seq_R1_processed.fastq
Approx 30% complete for K562_GRO-seq_R1_processed.fastq
Approx 35% complete for K562_GRO-seq_R1_processed.fastq
Approx 40% complete for K562_GRO-seq_R1_processed.fastq
Approx 45% complete for K562_GRO-seq_R1_processed.fastq
Approx 50% complete for K562_GRO-seq_R1_processed.fastq
Approx 55% complete for K562_GRO-seq_R1_processed.fastq
Approx 60% complete for K562_GRO-seq_R1_processed.fastq
Approx 65% complete for K562_GRO-seq_R1_processed.fastq
Approx 70% complete for K562_GRO-seq_R1_processed.fastq
Approx 75% complete for K562_GRO-seq_R1_processed.fastq
Approx 80% complete for K562_GRO-seq_R1_processed.fastq
Approx 85% complete for K562_GRO-seq_R1_processed.fastq
Approx 90% complete for K562_GRO-seq_R1_processed.fastq
Approx 95% complete for K562_GRO-seq_R1_processed.fastq
Analysis complete for K562_GRO-seq_R1_processed.fastq
#
# A fatal error has been detected by the Java Runtime Environment:
#
#  SIGSEGV (0xb) at pc=0x00007fec85b13dfc, pid=96637, tid=0x00007fec6ae36700
#
# JRE version: Java(TM) SE Runtime Environment (8.0_171-b11) (build 1.8.0_171-b11)
# Java VM: Java HotSpot(TM) 64-Bit Server VM (25.171-b11 mixed mode linux-amd64 compressed oops)
# Problematic frame:
# V  [libjvm.so+0x8e7dfc]  Monitor::wait(bool, long, bool)+0x36c
#
# Core dump written. Default location: /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/core or core.96637
#
# An error report file with more information is saved as:
# /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/hs_err_pid96637.log
#
# If you would like to submit a bug report, please visit:
#   http://bugreport.java.com/bugreport/crash.jsp
#
</pre>
Command completed. Elapsed time: 0:01:31. Running peak memory: 0.28GB.  
  PID: 96637;	Command: fastqc;	Return code: -6;	Memory used: 0.169GB

Subprocess returned nonzero result. Check above output for details
ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True
> `FastQC report r1`	fastqc/K562_GRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/processed_R1.flag` (96741)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.28GB.  
  PID: 96741;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (02-27 09:28:45) elapsed: 185.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt` (96742)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.28GB.  
  PID: 96742;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Adapter insertion distribution`	cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_GRO-seq_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	22	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:28:52) elapsed: 7.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/cutadapt/K562_GRO-seq_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.7189	PEPPRO	_RES_

### Prealignments (02-27 09:28:52) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:28:52) elapsed: 0.0 _TIME_


> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_GRO-seq -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/fastq/K562_GRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
29322367 reads; of these:
  29322367 (100.00%) were unpaired; of these:
    25012896 (85.30%) aligned 0 times
    4309471 (14.70%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
14.70% overall alignment rate

> `Aligned_reads_human_rDNA`	4309471.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	14.7	PEPPRO	_RES_

### Map to genome (02-27 09:31:30) elapsed: 157.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id K562_GRO-seq -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/tmpiiqfgiw1 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam` (97206,97207,97209)
<pre>
25012896 reads; of these:
  25012896 (100.00%) were unpaired; of these:
    2113226 (8.45%) aligned 0 times
    15215880 (60.83%) aligned exactly 1 time
    7683790 (30.72%) aligned >1 times
91.55% overall alignment rate
[bam_sort_core] merging from 6 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:10:33. Running peak memory: 3.594GB.  
  PID: 97207;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 97206;	Command: bowtie2;	Return code: 0;	Memory used: 3.594GB  
  PID: 97209;	Command: samtools;	Return code: 0;	Memory used: 0.888GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam` (98196)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 3.594GB.  
  PID: 98196;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	22899670	PEPPRO	_RES_

> `QC_filtered_reads`	4468262	PEPPRO	_RES_

> `Aligned_reads`	18431408	PEPPRO	_RES_

> `Alignment_rate`	62.86	PEPPRO	_RES_

> `Total_efficiency`	60.7	PEPPRO	_RES_

> `Read_depth`	2.65	PEPPRO	_RES_

### Compress all unmapped read files (02-27 09:52:31) elapsed: 1261.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/prealignments/K562_GRO-seq_human_rDNA_unmap.fq` (99183)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 3.594GB.  
  PID: 99183;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam` (99232)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.594GB.  
  PID: 99232;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	843399	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam` (99255)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.594GB.  
  PID: 99255;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_noMT.bam` (99270,99271,99272,99273)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.594GB.  
  PID: 99270;	Command: samtools;	Return code: 0;	Memory used: 0.001GB  
  PID: 99272;	Command: grep;	Return code: 0;	Memory used: 0.001GB  
  PID: 99271;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 99273;	Command: xargs;	Return code: 0;	Memory used: 0.043GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam` (99301)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 99301;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam` (99302)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.594GB.  
  PID: 99302;	Command: samtools;	Return code: 0;	Memory used: 0.016GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	50	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 09:54:56) elapsed: 146.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -c 8 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv` (99361)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/tmp_K562_GRO-seq_sort_ofjf347v'
Processing with 8 cores...
Discarding 91 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 104 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270511v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270374v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 3.594GB.  
  PID: 99361;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.365GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_bamQC.tsv`

> `NRF`	0.7	PEPPRO	_RES_

> `PBC1`	0.86	PEPPRO	_RES_

> `PBC2`	8.69	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_unmap.bam`  

> `samtools view -b -@ 8 -f 4  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_unmap.bam` (99605)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.594GB.  
  PID: 99605;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_temp.bam`

> `Unmapped_reads`	2113226	PEPPRO	_RES_

### Split BAM by strand (02-27 09:55:37) elapsed: 40.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam` (99638)
<pre>
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 3.594GB.  
  PID: 99638;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam` (99693)
<pre>
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 3.594GB.  
  PID: 99693;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 09:57:38) elapsed: 122.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (99747)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 99747;	Command: sed;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_plus_TssEnrichment.txt` (99749)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.594GB.  
  PID: 99749;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.37GB


> `TSS_coding_score`	23.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_minus_TssEnrichment.txt` (99773)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.594GB.  
  PID: 99773;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.438GB


> `TSS_non-coding_score`	8.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_minus_TssEnrichment.txt` (99797)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.594GB.  
  PID: 99797;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `TSS enrichment`	QC_hg38/K562_GRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_GRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt` (99817,99818,99819,99820)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 99817;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 99819;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 99818;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 99820;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt` (99822)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 99822;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-27 09:57:57) elapsed: 19.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_tss.bed` (99824,99825)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.594GB.  
  PID: 99824;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 99825;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (99829,99830)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 99829;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 99830;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSS_density.bed` (99832,99833,99834,99835)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.594GB.  
  PID: 99832;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB  
  PID: 99834;	Command: sort;	Return code: 0;	Memory used: 0.01GB  
  PID: 99833;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 99835;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_gene_body_density.bed` (99857,99858,99859)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 3.594GB.  
  PID: 99858;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 99857;	Command: bedtools;	Return code: 0;	Memory used: 0.052GB  
  PID: 99859;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed` (99888,99889,99890)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.594GB.  
  PID: 99888;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 99890;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 99889;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	14.58	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed` (99895)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.594GB.  
  PID: 99895;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `Pause index`	QC_hg38/K562_GRO-seq_pause_index.pdf	Pause index	QC_hg38/K562_GRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_pause_index.bed` (99915)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 99915;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 09:59:00) elapsed: 62.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam`
18431408 6290701

> `Plus_FRiP`	0.34	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam`
18431408 5983864

> `Minus_FRiP`	0.32	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_gene_sort.bed` (99965,99966)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.594GB.  
  PID: 99966;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 99965;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_gene_coverage.bed` (99969)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 3.594GB.  
  PID: 99969;	Command: bedtools;	Return code: 0;	Memory used: 0.075GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed.gz` (100228)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 100228;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed` (100229)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 100229;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 10:00:04) elapsed: 64.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/raw/hg38_annotations.bed` (100238)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.594GB.  
  PID: 100238;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR"` (100240)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 100240;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR_sort.bed` (100242,100243,100244,100245)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.594GB.  
  PID: 100242;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 100243;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 100245;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 100244;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_plus_coverage.bed` (100247)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.594GB.  
  PID: 100247;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_minus_coverage.bed` (100263)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.594GB.  
  PID: 100263;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR"` (100305)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 100305;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR_sort.bed` (100306,100307,100308,100309)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.594GB.  
  PID: 100306;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 100307;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 100309;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 100308;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_plus_coverage.bed` (100311)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.594GB.  
  PID: 100311;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_minus_coverage.bed` (100323)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.594GB.  
  PID: 100323;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer_sort.bed` (100334,100335,100336,100337)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.594GB.  
  PID: 100334;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 100335;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 100337;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB  
  PID: 100336;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_plus_coverage.bed` (100340)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.594GB.  
  PID: 100340;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_minus_coverage.bed` (100367)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.594GB.  
  PID: 100367;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Exon_sort.bed` (100379,100380,100381,100382)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.594GB.  
  PID: 100379;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 100380;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 100382;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 100381;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_plus_coverage.bed` (100390)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.594GB.  
  PID: 100390;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_minus_coverage.bed` (100402)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.594GB.  
  PID: 100402;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Intron_sort.bed` (100416,100417,100418,100419)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.594GB.  
  PID: 100416;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 100418;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 100417;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 100419;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_plus_coverage.bed` (100422)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.594GB.  
  PID: 100422;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_minus_coverage.bed` (100436)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.594GB.  
  PID: 100436;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_sort.bed` (100449,100450,100451,100452)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.594GB.  
  PID: 100449;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 100450;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 100452;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 100451;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_plus_coverage.bed` (100454)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.594GB.  
  PID: 100454;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_minus_coverage.bed` (100470)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.594GB.  
  PID: 100470;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region"` (100482)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 100482;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (100483,100484,100485,100486)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.594GB.  
  PID: 100483;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 100485;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 100484;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 100486;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (100489)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.594GB.  
  PID: 100489;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (100501)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.594GB.  
  PID: 100501;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


### Plot cFRiF/FRiF (02-27 10:03:08) elapsed: 184.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_GRO-seq -z 3099922541 -n 9058650 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (100523)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.594GB.  
  PID: 100523;	Command: Rscript;	Return code: 0;	Memory used: 0.526GB

> `cFRiF`	QC_hg38/K562_GRO-seq_cFRiF.pdf	cFRiF	QC_hg38/K562_GRO-seq_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_GRO-seq -z 3099922541 -n 9058650 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (100566)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.594GB.  
  PID: 100566;	Command: Rscript;	Return code: 0;	Memory used: 0.536GB

> `FRiF`	QC_hg38/K562_GRO-seq_FRiF.pdf	FRiF	QC_hg38/K562_GRO-seq_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 10:04:09) elapsed: 61.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_exons_sort.bed` (100599,100600)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.594GB.  
  PID: 100600;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB  
  PID: 100599;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_introns_sort.bed` (100607,100608,100609)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.594GB.  
  PID: 100608;	Command: bedtools;	Return code: 0;	Memory used: 0.078GB  
  PID: 100607;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 100609;	Command: bedtools;	Return code: 0;	Memory used: 0.092GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_coverage.bed` (100620)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.594GB.  
  PID: 100620;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_coverage.bed` (100641)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.594GB.  
  PID: 100641;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.431408)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_rpkm.bed` (100862,100863,100864)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.594GB.  
  PID: 100862;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 100864;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 100863;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.431408)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_rpkm.bed` (100867,100868,100869)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.594GB.  
  PID: 100867;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 100869;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 100868;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed` (100871,100872,100873)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 100871;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 100873;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 100872;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.33	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed --annotate` (100879)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.594GB.  
  PID: 100879;	Command: Rscript;	Return code: 0;	Memory used: 0.213GB

> `mRNA contamination`	QC_hg38/K562_GRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_GRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/QC_hg38/K562_GRO-seq_exon_intron_ratios.bed` (100902)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.594GB.  
  PID: 100902;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-27 10:05:20) elapsed: 72.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam` (100910)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.594GB.  
  PID: 100910;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_smooth_body_0-mer.bw -p 5 --variable-step` (100917)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_K562_GRO-seq_plus_cuttrace_39k1y94_'
Processing with 2 cores...
stdin is empty of data
stdin is empty of data
Discarding 100 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr14_KI270723v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 95 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270511v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270374v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 95 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_exact_body_0-mer.bw'
Merging 95 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:11:42. Running peak memory: 3.594GB.  
  PID: 100917;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.309GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam` (102825)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.594GB.  
  PID: 102825;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_smooth_body_0-mer.bw -p 5 --variable-step` (102831)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/aligned_hg38/K562_GRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_K562_GRO-seq_minus_cuttrace_0myyp153'
Processing with 2 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr11_KI270721v1_random', 'chr17_KI270730v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_exact_body_0-mer.bw'
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_GRO-seq/signal_hg38/K562_GRO-seq_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:11:00. Running peak memory: 3.594GB.  
  PID: 102831;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.544GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:03:44
*  Total elapsed time (all runs):  1:08:49
*         Peak memory (this run):  3.5943 GB
*        Pipeline completed time: 2020-02-27 10:28:17
