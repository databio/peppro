### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_40 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_40pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-aw29-25b
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/
*  Pipeline started at:   (06-11 17:13:57) elapsed: 7.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_40pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_40`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_40pctRNA.fastq.gz

> `File_mb`	5584.92	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:13:58) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_40pctRNA.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz` (388389)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 388389;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1.fastq` (388390)
<pre>
</pre>
Command completed. Elapsed time: 0:02:21. Running peak memory: 0.002GB.  
  PID: 388390;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	70000000	PEPPRO	_RES_

> `Fastq_reads`	70000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz']

### FASTQ processing:  (06-11 17:17:51) elapsed: 233.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_processed.fastq`  

> `(cutadapt -j 12 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt` (390042)
<pre>
</pre>
Command completed. Elapsed time: 0:02:05. Running peak memory: 4.637GB.  
  PID: 390042;	Command: cutadapt;	Return code: 0;	Memory used: 4.637GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_processed.fastq` (390564,390565)
<pre>
</pre>
Command completed. Elapsed time: 0:01:46. Running peak memory: 4.637GB.  
  PID: 390564;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 390565;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	45723195.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	1050230.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	1.5003	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastqc/K562_RNA-seq_40_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (391558)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.637GB.  
  PID: 391558;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	68949770	PEPPRO	_RES_

> `Trim_loss_rate`	1.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_processed.fastq` (391626)
<pre>
Started analysis of K562_RNA-seq_40_R1_processed.fastq
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:14) elapsed: 5843.0 _TIME_
Child process 391626 (fastqc) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/recover.lock.fastqc__K562_RNA-seq_40_R1_processed_fastqc.html
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/recover.lock.trimmed_fastqc

### Pipeline failed at:  (06-11 18:55:14) elapsed: 0.0 _TIME_

Total time: 1:41:25
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_40 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_40pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-32c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/
*  Pipeline started at:   (06-11 19:08:40) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_40pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_40`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_40pctRNA.fastq.gz

> `File_mb`	5584.92	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:08:40) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz']

### FASTQ processing:  (06-11 19:08:40) elapsed: 0.0 _TIME_


> `cutadapt --version`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_processed.fastq`  
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/lock.fastqc__K562_RNA-seq_40_R1_processed_fastqc.html
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastqc/K562_RNA-seq_40_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (82210)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 82210;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	68949770	PEPPRO	_RES_

> `Trim_loss_rate`	1.5	PEPPRO	_RES_
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/lock.trimmed_fastqc
Overwriting target...
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_processed.fastq` (82232)
<pre>
Started analysis of K562_RNA-seq_40_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_40_R1_processed.fastq
Analysis complete for K562_RNA-seq_40_R1_processed.fastq
#
# A fatal error has been detected by the Java Runtime Environment:
#
#  SIGSEGV (0xb) at pc=0x00007f3d15fe58a9, pid=82232, tid=0x00007f3cfaae7700
#
# JRE version: Java(TM) SE Runtime Environment (8.0_171-b11) (build 1.8.0_171-b11)
# Java VM: Java HotSpot(TM) 64-Bit Server VM (25.171-b11 mixed mode linux-amd64 compressed oops)
# Problematic frame:
# V  [libjvm.so+0x6e08a9]  jni_ReleasePrimitiveArrayCritical+0x69
#
# Core dump written. Default location: /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/core or core.82232
#
# An error report file with more information is saved as:
# /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/hs_err_pid82232.log
#
# If you would like to submit a bug report, please visit:
#   http://bugreport.java.com/bugreport/crash.jsp
#
</pre>
Command completed. Elapsed time: 0:03:56. Running peak memory: 0.243GB.  
  PID: 82232;	Command: fastqc;	Return code: -6;	Memory used: 0.243GB

Subprocess returned nonzero result. Check above output for details
ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True
> `FastQC report r1`	fastqc/K562_RNA-seq_40_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/processed_R1.flag` (82836)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.243GB.  
  PID: 82836;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-11 19:12:51) elapsed: 251.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt` (82837)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.243GB.  
  PID: 82837;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_40_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_40_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 19:12:57) elapsed: 6.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2309	PEPPRO	_RES_

### Prealignments (06-11 19:12:57) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:12:57) elapsed: 0.0 _TIME_


> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_40 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/prealignments/K562_RNA-seq_40_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
68949770 reads; of these:
  68949770 (100.00%) were unpaired; of these:
    64521780 (93.58%) aligned 0 times
    4427990 (6.42%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
6.42% overall alignment rate

> `Aligned_reads_human_rDNA`	4427990.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	6.42	PEPPRO	_RES_

### Map to genome (06-11 19:20:50) elapsed: 473.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam`  

> `bowtie2 -p 12 --very-sensitive --rg-id K562_RNA-seq_40 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/prealignments/K562_RNA-seq_40_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/tmph12anwx9 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam` (84828,84833,84834)
<pre>
64521780 reads; of these:
  64521780 (100.00%) were unpaired; of these:
    5066898 (7.85%) aligned 0 times
    41197278 (63.85%) aligned exactly 1 time
    18257604 (28.30%) aligned >1 times
92.15% overall alignment rate
[bam_sort_core] merging from 20 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:39:43. Running peak memory: 3.718GB.  
  PID: 84833;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 84828;	Command: bowtie2;	Return code: 0;	Memory used: 3.718GB  
  PID: 84834;	Command: samtools;	Return code: 0;	Memory used: 0.889GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam` (90348)
<pre>
</pre>
Command completed. Elapsed time: 0:02:39. Running peak memory: 3.718GB.  
  PID: 90348;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	59454882	PEPPRO	_RES_

> `QC_filtered_reads`	9022968	PEPPRO	_RES_

> `Aligned_reads`	50431914	PEPPRO	_RES_

> `Alignment_rate`	73.14	PEPPRO	_RES_

> `Total_efficiency`	72.05	PEPPRO	_RES_

> `Read_depth`	6.34	PEPPRO	_RES_

### Compress all unmapped read files (06-11 20:24:00) elapsed: 3789.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/prealignments/K562_RNA-seq_40_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/prealignments/K562_RNA-seq_40_human_rDNA_unmap.fq` (93635)
<pre>
</pre>
Command completed. Elapsed time: 0:02:24. Running peak memory: 3.718GB.  
  PID: 93635;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam` (93989)
<pre>
</pre>
Command completed. Elapsed time: 0:01:05. Running peak memory: 3.718GB.  
  PID: 93989;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	2452024	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam` (94164)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 3.718GB.  
  PID: 94164;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/chr_sizes.bed` (94209,94210,94211,94212)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.718GB.  
  PID: 94209;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 94211;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 94210;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 94212;	Command: grep;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_noMT.bam` (94215)
<pre>
</pre>
Command completed. Elapsed time: 0:01:04. Running peak memory: 3.718GB.  
  PID: 94215;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam` (94285)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.718GB.  
  PID: 94285;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam` (94287)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 3.718GB.  
  PID: 94287;	Command: samtools;	Return code: 0;	Memory used: 0.016GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 20:32:52) elapsed: 533.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv` (95006)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/tmp_K562_RNA-seq_40_sort_4vnlq6cu'
Processing with 12 cores...
Discarding 97 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 98 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:01:11. Running peak memory: 3.718GB.  
  PID: 95006;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.424GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv`

> `NRF`	0.75	PEPPRO	_RES_

> `PBC1`	0.9	PEPPRO	_RES_

> `PBC2`	13.91	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_unmap.bam`  

> `samtools view -b -@ 12 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_unmap.bam` (95120)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.718GB.  
  PID: 95120;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam`

> `Unmapped_reads`	5066898	PEPPRO	_RES_

### Split BAM by strand (06-11 20:34:31) elapsed: 98.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam` (95171)
<pre>
</pre>
Command completed. Elapsed time: 0:04:00. Running peak memory: 3.718GB.  
  PID: 95171;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam` (96009)
<pre>
</pre>
Command completed. Elapsed time: 0:03:52. Running peak memory: 3.718GB.  
  PID: 96009;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 20:42:23) elapsed: 472.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (96752)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.718GB.  
  PID: 96752;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_plus_TssEnrichment.txt` (96753)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.718GB.  
  PID: 96753;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.53GB


> `TSS_coding_score`	16.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_minus_TssEnrichment.txt` (96791)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.718GB.  
  PID: 96791;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.541GB


> `TSS_non-coding_score`	5.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_minus_TssEnrichment.txt` (96824)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.718GB.  
  PID: 96824;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_40_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_40_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt` (96846,96847,96848,96849)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.718GB.  
  PID: 96846;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 96848;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 96847;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 96849;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt` (96851)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.718GB.  
  PID: 96851;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-11 20:42:47) elapsed: 24.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_tss.bed` (96853,96854)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.718GB.  
  PID: 96853;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 96854;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_gene_body.bed` (96857,96858)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.718GB.  
  PID: 96857;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 96858;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSS_density.bed` (96861,96862,96863,96864)
<pre>
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 3.718GB.  
  PID: 96862;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 96864;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 96861;	Command: bedtools;	Return code: 0;	Memory used: 0.032GB  
  PID: 96863;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_gene_body_density.bed` (96930,96931,96932)
<pre>
</pre>
Command completed. Elapsed time: 0:01:46. Running peak memory: 3.718GB.  
  PID: 96930;	Command: bedtools;	Return code: 0;	Memory used: 0.176GB  
  PID: 96932;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 96931;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/tmpd3fxoegt` (97232,97233,97234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.718GB.  
  PID: 97232;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 97234;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 97233;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/tmpd3fxoegt | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 20:45:51) elapsed: 183.0 _TIME_

Total time: 1:37:11
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/tmpd3fxoegt | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_RNA-seq_40 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_40pctRNA.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-23a
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/
*  Pipeline started at:   (06-14 21:11:55) elapsed: 7.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_40pctRNA.fastq.gz']`
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
*        `sample_name`:  `K562_RNA-seq_40`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_40pctRNA.fastq.gz

> `File_mb`	5584.92	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:11:55) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz']

### FASTQ processing:  (06-14 21:11:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:11:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_40_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_40_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:11:55) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:11:55) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/prealignments/K562_RNA-seq_40_human_rDNA_unmap.fq

### Map to genome (06-14 21:11:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam`  

### Compress all unmapped read files (06-14 21:11:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/prealignments/K562_RNA-seq_40_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:11:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_unmap.bam`  

### Split BAM by strand (06-14 21:11:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam`  

### Calculate TSS enrichment (06-14 21:11:55) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_RNA-seq_40_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_40_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:11:55) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/tmp62k_foj0` (274384,274385,274386)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 274384;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 274386;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 274385;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/tmp62k_foj0 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0411612) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/tmp62k_foj0 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed` (274393)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 274393;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	8.63	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed` (274398)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.251GB.  
  PID: 274398;	Command: Rscript;	Return code: 0;	Memory used: 0.251GB

> `Pause index`	QC_hg38/K562_RNA-seq_40_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_40_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed` (274421)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.251GB.  
  PID: 274421;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:12:02) elapsed: 6.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam`
50431914 18750392

> `Plus_FRiP`	0.37	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam`
50431914 18315541

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_gene_sort.bed` (285553,285574)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.251GB.  
  PID: 285553;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 285574;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_gene_coverage.bed` (286228)
<pre>
</pre>
Command completed. Elapsed time: 0:01:36. Running peak memory: 0.251GB.  
  PID: 286228;	Command: bedtools;	Return code: 0;	Memory used: 0.126GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed.gz` (318793)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.251GB.  
  PID: 318793;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed` (318795)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.251GB.  
  PID: 318795;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:15:05) elapsed: 183.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed` (318803)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.251GB.  
  PID: 318803;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer_sort.bed` (318805,318806,318807,318808)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.251GB.  
  PID: 318805;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 318806;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 318808;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 318807;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_plus_coverage.bed` (318811)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 0.251GB.  
  PID: 318811;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_minus_coverage.bed` (318846)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 0.251GB.  
  PID: 318846;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_sort.bed` (318878,318879,318880,318881)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.251GB.  
  PID: 318878;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 318879;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 318881;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 318880;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_plus_coverage.bed` (318884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.251GB.  
  PID: 318884;	Command: bedtools;	Return code: 0;	Memory used: 0.088GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_minus_coverage.bed` (318914)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.251GB.  
  PID: 318914;	Command: bedtools;	Return code: 0;	Memory used: 0.066GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region"` (333241)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.251GB.  
  PID: 333241;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed` (333250,333251,333252,333253)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.251GB.  
  PID: 333250;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 333252;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 333251;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 333253;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_plus_coverage.bed` (334337)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 0.251GB.  
  PID: 334337;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_minus_coverage.bed` (351994)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.251GB.  
  PID: 351994;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR"` (363004)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.251GB.  
  PID: 363004;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR_sort.bed` (363006,363007,363008,363009)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.251GB.  
  PID: 363006;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 363007;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 363009;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 363008;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_plus_coverage.bed` (363011)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.251GB.  
  PID: 363011;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_minus_coverage.bed` (363089)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.251GB.  
  PID: 363089;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR"` (363143)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.251GB.  
  PID: 363143;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR_sort.bed` (363144,363145,363146,363147)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.251GB.  
  PID: 363144;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 363145;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 363147;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 363146;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_plus_coverage.bed` (363149)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 0.251GB.  
  PID: 363149;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_minus_coverage.bed` (363462)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.251GB.  
  PID: 363462;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon_sort.bed` (363492,363493,363494,363495)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.251GB.  
  PID: 363492;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 363493;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 363495;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 363494;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_plus_coverage.bed` (363500)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 0.251GB.  
  PID: 363500;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_minus_coverage.bed` (363535)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 0.251GB.  
  PID: 363535;	Command: bedtools;	Return code: 0;	Memory used: 0.043GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron_sort.bed` (374569,374570,374571,374572)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.251GB.  
  PID: 374569;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 374571;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 374570;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 374572;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_plus_coverage.bed` (374575)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 0.251GB.  
  PID: 374575;	Command: bedtools;	Return code: 0;	Memory used: 0.085GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_minus_coverage.bed` (389599)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 0.251GB.  
  PID: 389599;	Command: bedtools;	Return code: 0;	Memory used: 0.055GB


### Plot cFRiF/FRiF (06-14 21:23:51) elapsed: 526.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_40 -z 3099922541 -n 24418925 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_plus_coverage.bed` (407660)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.475GB.  
  PID: 407660;	Command: Rscript;	Return code: 0;	Memory used: 0.475GB

> `cFRiF`	QC_hg38/K562_RNA-seq_40_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_40_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_RNA-seq_40 -z 3099922541 -n 24418925 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_plus_coverage.bed` (407711)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 0.475GB.  
  PID: 407711;	Command: Rscript;	Return code: 0;	Memory used: 0.474GB

> `FRiF`	QC_hg38/K562_RNA-seq_40_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_40_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:24:54) elapsed: 63.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_exons_sort.bed` (407743,407744)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.475GB.  
  PID: 407743;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 407744;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_introns_sort.bed` (408815,408825,408826)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.475GB.  
  PID: 408815;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 408826;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 408825;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_coverage.bed` (411057)
<pre>
</pre>
Command completed. Elapsed time: 0:01:20. Running peak memory: 0.475GB.  
  PID: 411057;	Command: bedtools;	Return code: 0;	Memory used: 0.096GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_coverage.bed` (451940)
<pre>
</pre>
Command completed. Elapsed time: 0:01:18. Running peak memory: 0.475GB.  
  PID: 451940;	Command: bedtools;	Return code: 0;	Memory used: 0.091GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/50.431914)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_rpkm.bed` (452081,452082,452083)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.475GB.  
  PID: 452081;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 452083;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 452082;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/50.431914)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_rpkm.bed` (452085,452086,452087)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.475GB.  
  PID: 452085;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 452087;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 452086;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed` (452090,452091,452092)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.475GB.  
  PID: 452090;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 452092;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 452091;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	3.53	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed --annotate` (452098)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.475GB.  
  PID: 452098;	Command: Rscript;	Return code: 0;	Memory used: 0.317GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_40_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_40_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed` (452119)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.475GB.  
  PID: 452119;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-14 21:27:49) elapsed: 175.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam` (452128)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 0.475GB.  
  PID: 452128;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 50431914.0` (452147)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_40_plus_cuttrace_u697j3fp'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:57. Running peak memory: 2.598GB.  
  PID: 452147;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.598GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam` (39639)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 2.598GB.  
  PID: 39639;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 50431914.0` (39662)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_40_minus_cuttrace_9pp7snmr'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:52. Running peak memory: 2.6GB.  
  PID: 39662;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.6GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:32:36
*  Total elapsed time (all runs):  2:43:38
*         Peak memory (this run):  2.6003 GB
*        Pipeline completed time: 2020-06-14 21:44:23
