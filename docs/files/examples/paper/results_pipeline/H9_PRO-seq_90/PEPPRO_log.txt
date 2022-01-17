### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_90 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_90pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_90pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba26-35c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/
*  Pipeline started at:   (06-15 07:17:27) elapsed: 2.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_90pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_90pct_PE2.fastq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `12000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `H9_PRO-seq_90`
*              `scale`:  `True`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `PAIRED`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `8`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_90pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_90pct_PE2.fastq.gz

> `File_mb`	2146.02	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:28) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/H9_PRO-seq_90_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_90pct_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/H9_PRO-seq_90_R1.fastq.gz` (217612)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 217612;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/H9_PRO-seq_90_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/H9_PRO-seq_90_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_90pct_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/H9_PRO-seq_90_R2.fastq.gz` (217614)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 217614;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/H9_PRO-seq_90_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/H9_PRO-seq_90_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1.fastq` (217617)
<pre>
</pre>
Command completed. Elapsed time: 0:02:50. Running peak memory: 0.002GB.  
  PID: 217617;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/H9_PRO-seq_90_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2.fastq` (350867)
<pre>
</pre>
Command completed. Elapsed time: 0:02:27. Running peak memory: 0.002GB.  
  PID: 350867;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	104146252	PEPPRO	_RES_

> `Fastq_reads`	104146252	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/H9_PRO-seq_90_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/H9_PRO-seq_90_R2.fastq.gz']

### FASTQ processing:  (06-15 07:26:06) elapsed: 518.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90_R1_cutadapt.txt` (332572)
<pre>
</pre>
Command completed. Elapsed time: 0:02:48. Running peak memory: 3.094GB.  
  PID: 332572;	Command: cutadapt;	Return code: 0;	Memory used: 3.094GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq` (332770,332771)
<pre>
</pre>
Command completed. Elapsed time: 0:00:52. Running peak memory: 3.094GB.  
  PID: 332770;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 332771;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	25498419	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq` (332870)
<pre>
Started analysis of H9_PRO-seq_90_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_90_R1_processed.fastq
Analysis complete for H9_PRO-seq_90_R1_processed.fastq
#
# A fatal error has been detected by the Java Runtime Environment:
#
#  SIGSEGV (0xb) at pc=0x00007fba053fb8a9, pid=332870, tid=0x00007fb9e9eec700
#
# JRE version: Java(TM) SE Runtime Environment (8.0_171-b11) (build 1.8.0_171-b11)
# Java VM: Java HotSpot(TM) 64-Bit Server VM (25.171-b11 mixed mode linux-amd64 compressed oops)
# Problematic frame:
# V  [libjvm.so+0x6e08a9]  jni_ReleasePrimitiveArrayCritical+0x69
#
# Core dump written. Default location: /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/core or core.332870
#
# An error report file with more information is saved as:
# /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/hs_err_pid332870.log
#
# If you would like to submit a bug report, please visit:
#   http://bugreport.java.com/bugreport/crash.jsp
#
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 3.094GB.  
  PID: 332870;	Command: fastqc;	Return code: -6;	Memory used: 0.178GB

Subprocess returned nonzero result. Check above output for details
ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True
> `FastQC report r1`	fastqc/H9_PRO-seq_90_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_noadap.fastq` (333200)
<pre>
[INFO][0m 2562794 duplicated records removed
</pre>
Command completed. Elapsed time: 0:01:02. Running peak memory: 3.094GB.  
  PID: 333200;	Command: seqkit;	Return code: 0;	Memory used: 2.048GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_trimmed.fastq` (333284,333285)
<pre>
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 3.094GB.  
  PID: 333284;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 333285;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	41262698.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	25749458.0	PEPPRO	_RES_

> `Duplicate_reads`	2562794.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.4487	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/processed_R1.flag` (333542)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.094GB.  
  PID: 333542;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90_R2_cutadapt.txt` (333544)
<pre>
</pre>
Command completed. Elapsed time: 0:02:09. Running peak memory: 3.381GB.  
  PID: 333544;	Command: cutadapt;	Return code: 0;	Memory used: 3.381GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed.fastq` (333910,333911)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 3.381GB.  
  PID: 333910;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 333911;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	50996838	PEPPRO	_RES_

> `Trim_loss_rate`	51.03	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq` (334000)
<pre>
Started analysis of H9_PRO-seq_90_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_90_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_90_R1_processed.fastq
Analysis complete for H9_PRO-seq_90_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 3.381GB.  
  PID: 334000;	Command: fastqc;	Return code: 0;	Memory used: 0.179GB

> `FastQC report r1`	fastqc/H9_PRO-seq_90_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed.fastq` (334070)
<pre>
Started analysis of H9_PRO-seq_90_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_90_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_90_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_90_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:01:05. Running peak memory: 3.381GB.  
  PID: 334070;	Command: fastqc;	Return code: 0;	Memory used: 0.166GB

> `FastQC report r2`	fastqc/H9_PRO-seq_90_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90.histogram`  

> `fastq_pair -t 93731626 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_noadap.fastq` (334203)
<pre>
Left paired: 26141450		Right paired: 26141450
Left single: 182218		Right single: 2409759
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:03:26. Running peak memory: 5.91GB.  
  PID: 334203;	Command: fastq_pair;	Return code: 0;	Memory used: 5.91GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_noadap.fastq.paired.fq -o H9_PRO-seq_90 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt` (334807)
<pre>
</pre>
Command completed. Elapsed time: 0:01:24. Running peak memory: 5.91GB.  
  PID: 334807;	Command: flash;	Return code: 0;	Memory used: 0.126GB


### Plot adapter insertion distribution (06-15 07:44:33) elapsed: 1107.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt -u 8` (335026)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.91GB.  
  PID: 335026;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_90_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_90_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:44:38) elapsed: 5.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/cutadapt/H9_PRO-seq_90.hist`

> `Degradation_ratio`	1.0318	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed_dups.fastq` (335056)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 5.91GB.  
  PID: 335056;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/processed_R2.flag` (335273)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 335273;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/repaired.flag`  

> `fastq_pair -t 93731626 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed.fastq` (335274)
<pre>
Left paired: 25297086		Right paired: 25297086
Left single: 201333		Right single: 2474476
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:55. Running peak memory: 5.91GB.  
  PID: 335274;	Command: fastq_pair;	Return code: 0;	Memory used: 5.758GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq` (335620)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 335620;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed.fastq` (335621)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 335621;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/repaired.flag` (335623)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 335623;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/dups_repaired.flag`  

> `fastq_pair -t 93731626 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed_dups.fastq` (335624)
<pre>
Left paired: 22929753		Right paired: 22929753
Left single: 147931		Right single: 4841809
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:19. Running peak memory: 5.91GB.  
  PID: 335624;	Command: fastq_pair;	Return code: 0;	Memory used: 5.254GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_trimmed.fastq` (335998)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 335998;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed_dups.fastq` (336000)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 336000;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/dups_repaired.flag` (336001)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 336001;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:50:26) elapsed: 348.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:50:27) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/human_rDNA_bt2` (336003)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 336003;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_R2.fq` (336004)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_90 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
25297086 reads; of these:
  25297086 (100.00%) were unpaired; of these:
    22674298 (89.63%) aligned 0 times
    2622788 (10.37%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.37% overall alignment rate

> `Aligned_reads_human_rDNA`	5245576.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.29	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:54:33) elapsed: 246.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/human_rDNA_dups_bt2` (336273)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 336273;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_dups_R2.fq` (336274)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_90 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/fastq/H9_PRO-seq_90_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2622788 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:58:04) elapsed: 211.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_90 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/tmp6n7q80nc -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp.bam` (336675,336676,336677)
<pre>
2051450 reads skipped
0 reads lost
22674298 reads; of these:
  22674298 (100.00%) were paired; of these:
    9399583 (41.45%) aligned concordantly 0 times
    10886810 (48.01%) aligned concordantly exactly 1 time
    2387905 (10.53%) aligned concordantly >1 times
    ----
    9399583 pairs aligned concordantly 0 times; of these:
      2091316 (22.25%) aligned discordantly 1 time
    ----
    7308267 pairs aligned 0 times concordantly or discordantly; of these:
      14616534 mates make up the pairs; of these:
        5998097 (41.04%) aligned 0 times
        3149397 (21.55%) aligned exactly 1 time
        5469040 (37.42%) aligned >1 times
86.77% overall alignment rate
[bam_sort_core] merging from 12 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:46:22. Running peak memory: 5.91GB.  
  PID: 336676;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 336675;	Command: bowtie2;	Return code: 0;	Memory used: 3.788GB  
  PID: 336677;	Command: samtools;	Return code: 0;	Memory used: 0.913GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort.bam` (341319)
<pre>
</pre>
Command completed. Elapsed time: 0:01:23. Running peak memory: 5.91GB.  
  PID: 341319;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	39350499	PEPPRO	_RES_

> `QC_filtered_reads`	23000545	PEPPRO	_RES_

> `Aligned_reads`	16349954.5	PEPPRO	_RES_

> `Alignment_rate`	32.06	PEPPRO	_RES_

> `Total_efficiency`	15.7	PEPPRO	_RES_

> `Read_depth`	2.95	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_90 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/tmp6n7q80nc -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp_dups.bam` (342692,342697,342698)
<pre>
20878303 reads; of these:
  20878303 (100.00%) were paired; of these:
    8468903 (40.56%) aligned concordantly 0 times
    10204244 (48.87%) aligned concordantly exactly 1 time
    2205156 (10.56%) aligned concordantly >1 times
    ----
    8468903 pairs aligned concordantly 0 times; of these:
      1966476 (23.22%) aligned discordantly 1 time
    ----
    6502427 pairs aligned 0 times concordantly or discordantly; of these:
      13004854 mates make up the pairs; of these:
        5435246 (41.79%) aligned 0 times
        2954556 (22.72%) aligned exactly 1 time
        4615052 (35.49%) aligned >1 times
86.98% overall alignment rate
[bam_sort_core] merging from 11 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:38:37. Running peak memory: 5.91GB.  
  PID: 342692;	Command: bowtie2;	Return code: 0;	Memory used: 3.762GB  
  PID: 342697;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 342698;	Command: samtools;	Return code: 0;	Memory used: 0.977GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort_dups.bam` (346945)
<pre>
</pre>
Command completed. Elapsed time: 0:01:19. Running peak memory: 5.91GB.  
  PID: 346945;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 09:38:22) elapsed: 6018.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_R2.fq` (347073)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 5.91GB.  
  PID: 347073;	Command: pigz;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/prealignments/H9_PRO-seq_90_human_rDNA_unmap_R1.fq` (347114)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 5.91GB.  
  PID: 347114;	Command: pigz;	Return code: 0;	Memory used: 0.013GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp.bam` (347156)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 5.91GB.  
  PID: 347156;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	756629	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort.bam` (347236)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 5.91GB.  
  PID: 347236;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/chr_sizes.bed` (347282,347283,347284,347285)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 347283;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 347285;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 347282;	Command: samtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 347284;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_noMT.bam` (347287)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 5.91GB.  
  PID: 347287;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort.bam` (347574)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 347574;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort.bam` (347575)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 5.91GB.  
  PID: 347575;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


### Split BAM file (06-15 09:40:56) elapsed: 154.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam` (347602,347603)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:17. Running peak memory: 5.91GB.  
  PID: 347603;	Command: samtools;	Return code: 0;	Memory used: 4.627GB  
  PID: 347602;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE2.bam` (347924,347925)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:52. Running peak memory: 5.91GB.  
  PID: 347924;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 347925;	Command: samtools;	Return code: 0;	Memory used: 3.745GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp_dups.bam` (348964)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 5.91GB.  
  PID: 348964;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_dups_PE1.bam` (348998,348999)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:09. Running peak memory: 5.91GB.  
  PID: 348998;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 348999;	Command: samtools;	Return code: 0;	Memory used: 4.465GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_dups_PE2.bam` (349202,349203)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:52. Running peak memory: 5.91GB.  
  PID: 349202;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 349203;	Command: samtools;	Return code: 0;	Memory used: 3.58GB


### Calculate library complexity (06-15 09:50:27) elapsed: 571.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_dups_PE1.bam` (349728)
<pre>
BAM_INPUT
TOTAL READS     = 16739307
COUNTS_SUM      = 16739307
DISTINCT READS  = 1.35944e+07
DISTINCT COUNTS = 302
MAX COUNT       = 25272
COUNTS OF 1     = 1.21094e+07
OBSERVED COUNTS (25273)
1	12109447
2	1009536
3	235059
4	92617
5	47577
6	27424
7	17453
8	11457
9	8401
10	6126
11	4555
12	3543
13	2867
14	2217
15	1761
16	1480
17	1206
18	1066
19	933
20	821
21	720
22	595
23	543
24	462
25	433
26	358
27	352
28	312
29	298
30	244
31	251
32	210
33	209
34	181
35	151
36	142
37	135
38	121
39	138
40	124
41	119
42	133
43	91
44	103
45	87
46	88
47	88
48	81
49	63
50	65
51	75
52	68
53	69
54	59
55	60
56	52
57	50
58	46
59	36
60	28
61	48
62	39
63	46
64	46
65	21
66	31
67	35
68	26
69	29
70	28
71	30
72	19
73	28
74	30
75	25
76	18
77	27
78	21
79	21
80	17
81	22
82	13
83	16
84	24
85	20
86	12
87	20
88	9
89	15
90	22
91	15
92	18
93	12
94	8
95	4
96	11
97	13
98	14
99	12
100	13
101	7
102	10
103	11
104	7
105	6
106	9
107	9
108	9
109	7
110	5
111	7
112	11
113	10
114	9
115	4
116	6
117	7
118	5
119	6
121	10
122	9
123	5
124	5
125	10
126	3
127	5
128	3
129	8
130	2
131	8
132	6
133	6
134	4
135	2
136	3
137	2
138	2
139	4
140	6
141	5
142	6
143	4
144	6
145	5
146	7
147	1
148	3
149	2
150	4
151	1
152	4
153	3
154	4
156	5
157	2
158	4
159	3
161	1
162	3
163	2
164	3
165	5
166	1
167	1
168	1
170	1
171	2
172	1
173	2
174	3
176	3
178	5
179	4
180	1
181	3
182	4
183	2
184	2
185	3
186	3
187	1
188	2
190	3
192	2
193	2
194	4
195	2
196	2
198	3
199	1
200	1
201	4
202	1
204	1
205	1
206	2
207	2
208	3
209	2
210	3
212	1
213	3
214	1
215	1
216	1
219	1
223	2
226	2
231	1
234	2
236	4
237	1
238	1
241	1
242	1
243	2
244	1
245	1
247	2
248	2
250	1
251	1
254	2
262	1
267	1
271	1
272	1
274	1
276	1
280	1
281	1
282	1
283	1
285	1
286	2
288	1
289	2
291	1
292	1
294	2
295	2
296	1
300	1
302	1
310	1
311	2
317	1
325	1
326	1
332	1
333	1
334	1
341	1
344	1
348	1
355	1
357	2
362	1
367	1
389	1
392	1
398	1
402	1
411	1
445	1
448	1
453	1
459	1
461	2
466	1
490	1
493	1
495	1
529	1
623	1
637	1
643	1
646	2
648	1
649	1
772	1
800	1
828	1
878	1
910	1
1059	1
1140	1
1197	1
1468	1
1699	1
1712	1
1716	1
2073	1
2212	1
2338	1
3114	1
4028	1
8278	1
9628	1
12271	1
21870	1
25272	1

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
sample size: 14000000
sample size: 15000000
sample size: 16000000
</pre>
Command completed. Elapsed time: 0:01:30. Running peak memory: 5.91GB.  
  PID: 349728;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_dups_PE1.bam` (350400)
<pre>
BAM_INPUT
TOTAL READS     = 16739307
DISTINCT READS  = 1.35944e+07
DISTINCT COUNTS = 302
MAX COUNT       = 25272
COUNTS OF 1     = 1.21094e+07
MAX TERMS       = 100
OBSERVED COUNTS (25273)
1	12109447
2	1009536
3	235059
4	92617
5	47577
6	27424
7	17453
8	11457
9	8401
10	6126
11	4555
12	3543
13	2867
14	2217
15	1761
16	1480
17	1206
18	1066
19	933
20	821
21	720
22	595
23	543
24	462
25	433
26	358
27	352
28	312
29	298
30	244
31	251
32	210
33	209
34	181
35	151
36	142
37	135
38	121
39	138
40	124
41	119
42	133
43	91
44	103
45	87
46	88
47	88
48	81
49	63
50	65
51	75
52	68
53	69
54	59
55	60
56	52
57	50
58	46
59	36
60	28
61	48
62	39
63	46
64	46
65	21
66	31
67	35
68	26
69	29
70	28
71	30
72	19
73	28
74	30
75	25
76	18
77	27
78	21
79	21
80	17
81	22
82	13
83	16
84	24
85	20
86	12
87	20
88	9
89	15
90	22
91	15
92	18
93	12
94	8
95	4
96	11
97	13
98	14
99	12
100	13
101	7
102	10
103	11
104	7
105	6
106	9
107	9
108	9
109	7
110	5
111	7
112	11
113	10
114	9
115	4
116	6
117	7
118	5
119	6
121	10
122	9
123	5
124	5
125	10
126	3
127	5
128	3
129	8
130	2
131	8
132	6
133	6
134	4
135	2
136	3
137	2
138	2
139	4
140	6
141	5
142	6
143	4
144	6
145	5
146	7
147	1
148	3
149	2
150	4
151	1
152	4
153	3
154	4
156	5
157	2
158	4
159	3
161	1
162	3
163	2
164	3
165	5
166	1
167	1
168	1
170	1
171	2
172	1
173	2
174	3
176	3
178	5
179	4
180	1
181	3
182	4
183	2
184	2
185	3
186	3
187	1
188	2
190	3
192	2
193	2
194	4
195	2
196	2
198	3
199	1
200	1
201	4
202	1
204	1
205	1
206	2
207	2
208	3
209	2
210	3
212	1
213	3
214	1
215	1
216	1
219	1
223	2
226	2
231	1
234	2
236	4
237	1
238	1
241	1
242	1
243	2
244	1
245	1
247	2
248	2
250	1
251	1
254	2
262	1
267	1
271	1
272	1
274	1
276	1
280	1
281	1
282	1
283	1
285	1
286	2
288	1
289	2
291	1
292	1
294	2
295	2
296	1
300	1
302	1
310	1
311	2
317	1
325	1
326	1
332	1
333	1
334	1
341	1
344	1
348	1
355	1
357	2
362	1
367	1
389	1
392	1
398	1
402	1
411	1
445	1
448	1
453	1
459	1
461	2
466	1
490	1
493	1
495	1
529	1
623	1
637	1
643	1
646	2
648	1
649	1
772	1
800	1
828	1
878	1
910	1
1059	1
1140	1
1197	1
1468	1
1699	1
1712	1
1716	1
2073	1
2212	1
2338	1
3114	1
4028	1
8278	1
9628	1
12271	1
21870	1
25272	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
.....__..__......_.._................._..........__......_...._....._._...___....._..._.._........_................_.....
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:45. Running peak memory: 5.91GB.  
  PID: 350400;	Command: preseq;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_counts.txt` (350557)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 5.91GB.  
  PID: 350557;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_plot` (350611)
<pre>
Processing H9_PRO-seq_90
INFO: Found real counts for H9_PRO-seq_90 - Total (M): 17.534903 Unique (M): 16.739307

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.91GB.  
  PID: 350611;	Command: Rscript;	Return code: 0;	Memory used: 0.32GB

> `Library complexity`	QC_hg38/H9_PRO-seq_90_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_90_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8542	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 09:54:11) elapsed: 224.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam` (350630)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 5.91GB.  
  PID: 350630;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_bamQC.tsv` (350642)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/tmp_H9_PRO-seq_90_PE1_eyt9l283'
Processing with 12 cores...
Discarding 100 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270755v1']
Keeping 95 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 5.91GB.  
  PID: 350642;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.599GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	8767451.5	PEPPRO	_RES_

> `PBC2`	8767451.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_unmap.bam` (350688)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 5.91GB.  
  PID: 350688;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_temp.bam`

> `Unmapped_reads`	5998097	PEPPRO	_RES_

### Split BAM by strand (06-15 09:54:53) elapsed: 42.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam` (350731)
<pre>
</pre>
Command completed. Elapsed time: 0:01:02. Running peak memory: 5.91GB.  
  PID: 350731;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam` (351023)
<pre>
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 5.91GB.  
  PID: 351023;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 09:56:54) elapsed: 122.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (351103)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351103;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_plus_TssEnrichment.txt` (351105)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 5.91GB.  
  PID: 351105;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.971GB


> `TSS_coding_score`	33.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_minus_TssEnrichment.txt` (351136)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.91GB.  
  PID: 351136;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.949GB


> `TSS_non-coding_score`	12.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_minus_TssEnrichment.txt` (351169)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 5.91GB.  
  PID: 351169;	Command: Rscript;	Return code: 0;	Memory used: 0.319GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_90_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_90_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt` (351190,351191,351192,351193)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351190;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 351192;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 351191;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 351193;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_keep.txt` (351195)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351195;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 09:57:12) elapsed: 17.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_ensembl_tss.bed` (351197,351198)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 5.91GB.  
  PID: 351197;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 351198;	Command: bedtools;	Return code: 0;	Memory used: 0.095GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_ensembl_gene_body.bed` (351202,351203)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351202;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 351203;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_TSS_density.bed` (351205,351206,351207,351208)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 5.91GB.  
  PID: 351205;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 351207;	Command: sort;	Return code: 0;	Memory used: 0.012GB  
  PID: 351206;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 351208;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_gene_body_density.bed` (351228,351229,351230)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 5.91GB.  
  PID: 351229;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 351228;	Command: bedtools;	Return code: 0;	Memory used: 0.09GB  
  PID: 351230;	Command: sort;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/tmpf1u59rxb` (351277,351278,351279)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351277;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 351279;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 351278;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/tmpf1u59rxb | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0171122) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/tmpf1u59rxb > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_pause_index.bed` (351285)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351285;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	33.96	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_pause_index.bed` (351290)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.91GB.  
  PID: 351290;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/H9_PRO-seq_90_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_90_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_pause_index.bed` (351311)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351311;	Command: pigz;	Return code: 0;	Memory used: 0.004GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 09:58:07) elapsed: 55.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam`
16349954.5 6344819

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam`
16349954.5 5914666

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_gene_sort.bed` (351355,351356)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 351355;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 351356;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_gene_coverage.bed` (351359)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 5.91GB.  
  PID: 351359;	Command: bedtools;	Return code: 0;	Memory used: 0.092GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/hg38_annotations.bed.gz` (351380)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351380;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/hg38_annotations.bed` (351381)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351381;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 09:59:04) elapsed: 57.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/raw/hg38_annotations.bed` (351390)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 351390;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Enhancer_sort.bed` (351392,351393,351394,351395)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 351392;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 351393;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 351395;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 351394;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Enhancer_plus_coverage.bed` (351398)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 5.91GB.  
  PID: 351398;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Enhancer_minus_coverage.bed` (351408)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 5.91GB.  
  PID: 351408;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_sort.bed` (351418,351419,351420,351421)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351418;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 351419;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 351421;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 351420;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_plus_coverage.bed` (351424)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 5.91GB.  
  PID: 351424;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_minus_coverage.bed` (351436)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 5.91GB.  
  PID: 351436;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_Flanking_Region"` (351447)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351447;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_Flanking_Region_sort.bed` (351448,351449,351450,351451)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 351448;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 351450;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 351449;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 351451;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_Flanking_Region_plus_coverage.bed` (351454)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 5.91GB.  
  PID: 351454;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_Flanking_Region_minus_coverage.bed` (351472)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 5.91GB.  
  PID: 351472;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/5_UTR"` (351722)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351722;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/5_UTR_sort.bed` (351723,351724,351725,351726)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 351723;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 351724;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 351726;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 351725;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_5_UTR_plus_coverage.bed` (351728)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 5.91GB.  
  PID: 351728;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_5_UTR_minus_coverage.bed` (351739)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 5.91GB.  
  PID: 351739;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/3_UTR"` (351749)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 351749;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/3_UTR_sort.bed` (351750,351751,351752,351753)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 351750;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 351751;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 351753;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 351752;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_3_UTR_plus_coverage.bed` (351756)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 5.91GB.  
  PID: 351756;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_3_UTR_minus_coverage.bed` (351766)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 5.91GB.  
  PID: 351766;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Exon_sort.bed` (351778,351779,351780,351781)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 5.91GB.  
  PID: 351778;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 351779;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 351781;	Command: bedtools;	Return code: 0;	Memory used: 0.172GB  
  PID: 351780;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Exon_plus_coverage.bed` (351785)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 5.91GB.  
  PID: 351785;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Exon_minus_coverage.bed` (351808)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 5.91GB.  
  PID: 351808;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Intron_sort.bed` (351820,351821,351822,351823)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 351820;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 351822;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 351821;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 351823;	Command: bedtools;	Return code: 0;	Memory used: 0.082GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Intron_plus_coverage.bed` (351829)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.91GB.  
  PID: 351829;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Intron_minus_coverage.bed` (351844)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 5.91GB.  
  PID: 351844;	Command: bedtools;	Return code: 0;	Memory used: 0.031GB


### Plot cFRiF/FRiF (06-15 10:01:48) elapsed: 165.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_90 -z 3099922541 -n 8950537 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Intron_plus_coverage.bed` (351882)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 5.91GB.  
  PID: 351882;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `cFRiF`	QC_hg38/H9_PRO-seq_90_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_90_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_90 -z 3099922541 -n 8950537 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_Intron_plus_coverage.bed` (351973)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 5.91GB.  
  PID: 351973;	Command: Rscript;	Return code: 0;	Memory used: 0.458GB

> `FRiF`	QC_hg38/H9_PRO-seq_90_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_90_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 10:02:48) elapsed: 60.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_exons_sort.bed` (352021,352022)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.91GB.  
  PID: 352022;	Command: bedtools;	Return code: 0;	Memory used: 0.1GB  
  PID: 352021;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_introns_sort.bed` (352028,352029,352030)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.91GB.  
  PID: 352028;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 352030;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 352029;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exons_coverage.bed` (352036)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 5.91GB.  
  PID: 352036;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_introns_coverage.bed` (352054)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 5.91GB.  
  PID: 352054;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/16.3499545)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exons_rpkm.bed` (352079,352080,352081)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 352079;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 352081;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 352080;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/16.3499545)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_introns_rpkm.bed` (352083,352084,352085)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.91GB.  
  PID: 352083;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 352085;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 352084;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exon_intron_ratios.bed` (352088,352089,352090)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 352088;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 352090;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 352089;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exon_intron_ratios.bed --annotate` (352096)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.91GB.  
  PID: 352096;	Command: Rscript;	Return code: 0;	Memory used: 0.317GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_90_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_90_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/QC_hg38/H9_PRO-seq_90_exon_intron_ratios.bed` (352116)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.91GB.  
  PID: 352116;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 10:03:48) elapsed: 60.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam` (352125)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 5.91GB.  
  PID: 352125;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 16349954.5` (352131)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_90_plus_cuttrace_rfo82le1'
Processing with 4 cores...
stdin is empty of data
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_plus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:01. Running peak memory: 5.91GB.  
  PID: 352131;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.6GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam` (353975)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 5.91GB.  
  PID: 353975;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 16349954.5` (353985)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/aligned_hg38/H9_PRO-seq_90_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_90_minus_cuttrace_huwmlmjq'
Processing with 4 cores...
stdin is empty of data
Discarding 113 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr14_KI270723v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270755v1']
Keeping 82 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_minus_exact_body_0-mer.bw'
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_90/signal_hg38/H9_PRO-seq_90_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:13. Running peak memory: 5.91GB.  
  PID: 353985;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.361GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:00:51
*  Total elapsed time (all runs):  5:15:08
*         Peak memory (this run):  5.9105 GB
*        Pipeline completed time: 2020-06-15 10:18:16
