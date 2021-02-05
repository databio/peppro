### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_70 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_70pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_70pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-33c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/
*  Pipeline started at:   (06-15 07:17:18) elapsed: 1.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_70pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_70pct_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_70`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_70pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_70pct_PE2.fastq.gz

> `File_mb`	1677.91	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:20) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/H9_PRO-seq_70_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_70pct_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/H9_PRO-seq_70_R1.fastq.gz` (420358)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 420358;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/H9_PRO-seq_70_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/H9_PRO-seq_70_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_70pct_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/H9_PRO-seq_70_R2.fastq.gz` (420361)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 420361;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/H9_PRO-seq_70_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/H9_PRO-seq_70_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1.fastq` (420364)
<pre>
</pre>
Command completed. Elapsed time: 0:02:04. Running peak memory: 0.002GB.  
  PID: 420364;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/H9_PRO-seq_70_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2.fastq` (420497)
<pre>
</pre>
Command completed. Elapsed time: 0:01:35. Running peak memory: 0.002GB.  
  PID: 420497;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	81008246	PEPPRO	_RES_

> `Fastq_reads`	81008246	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/H9_PRO-seq_70_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/H9_PRO-seq_70_R2.fastq.gz']

### FASTQ processing:  (06-15 07:22:29) elapsed: 310.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70_R1_cutadapt.txt` (421060)
<pre>
</pre>
Command completed. Elapsed time: 0:01:26. Running peak memory: 3.108GB.  
  PID: 421060;	Command: cutadapt;	Return code: 0;	Memory used: 3.108GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq` (421227,421228)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.108GB.  
  PID: 421227;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 421228;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	19833681	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq` (421286)
<pre>
Started analysis of H9_PRO-seq_70_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_70_R1_processed.fastq
Analysis complete for H9_PRO-seq_70_R1_processed.fastq
#
# A fatal error has been detected by the Java Runtime Environment:
#
#  SIGSEGV (0xb) at pc=0x00007f91aa8a2bd4, pid=421286, tid=0x00007f91798d3700
#
# JRE version: Java(TM) SE Runtime Environment (8.0_171-b11) (build 1.8.0_171-b11)
# Java VM: Java HotSpot(TM) 64-Bit Server VM (25.171-b11 mixed mode linux-amd64 compressed oops)
# Problematic frame:
# V  [libjvm.so+0x418bd4]  ciEnv::get_field_by_index(ciInstanceKlass*, int)+0x94
#
# Core dump written. Default location: /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/core or core.421286
#
# An error report file with more information is saved as:
# /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/hs_err_pid421286.log
[thread 140262776989440 also had an error]
# [ timer expired, abort... ]
</pre>
Command completed. Elapsed time: 0:02:42. Running peak memory: 3.108GB.  
  PID: 421286;	Command: fastqc;	Return code: -6;	Memory used: 0.216GB

Subprocess returned nonzero result. Check above output for details
ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True
> `FastQC report r1`	fastqc/H9_PRO-seq_70_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_noadap.fastq` (421765)
<pre>
[INFO][0m 1700774 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 3.108GB.  
  PID: 421765;	Command: seqkit;	Return code: 0;	Memory used: 1.035GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_trimmed.fastq` (421810,421811)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.108GB.  
  PID: 421810;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 421811;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	32096850.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	20028832.0	PEPPRO	_RES_

> `Duplicate_reads`	1700774.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.4489	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/processed_R1.flag` (422058)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.108GB.  
  PID: 422058;	Command: touch;	Return code: 0;	Memory used: 0.002GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70_R2_cutadapt.txt` (422060)
<pre>
</pre>
Command completed. Elapsed time: 0:01:11. Running peak memory: 3.416GB.  
  PID: 422060;	Command: cutadapt;	Return code: 0;	Memory used: 3.416GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed.fastq` (422176,422177)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.416GB.  
  PID: 422176;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 422177;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	39667362	PEPPRO	_RES_

> `Trim_loss_rate`	51.03	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq` (422468)
<pre>
Started analysis of H9_PRO-seq_70_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_70_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_70_R1_processed.fastq
Analysis complete for H9_PRO-seq_70_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.416GB.  
  PID: 422468;	Command: fastqc;	Return code: 0;	Memory used: 0.166GB

> `FastQC report r1`	fastqc/H9_PRO-seq_70_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed.fastq` (422521)
<pre>
Started analysis of H9_PRO-seq_70_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_70_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_70_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_70_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 3.416GB.  
  PID: 422521;	Command: fastqc;	Return code: 0;	Memory used: 0.17GB

> `FastQC report r2`	fastqc/H9_PRO-seq_70_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70.histogram`  

> `fastq_pair -t 72907421 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_noadap.fastq` (422735)
<pre>
Left paired: 20333657		Right paired: 20333657
Left single: 141634		Right single: 1875182
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:12. Running peak memory: 4.435GB.  
  PID: 422735;	Command: fastq_pair;	Return code: 0;	Memory used: 4.435GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_noadap.fastq.paired.fq -o H9_PRO-seq_70 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt` (422905)
<pre>
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 4.435GB.  
  PID: 422905;	Command: flash;	Return code: 0;	Memory used: 0.116GB


### Plot adapter insertion distribution (06-15 07:34:46) elapsed: 737.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt -u 8` (423091)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.435GB.  
  PID: 423091;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_70_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_70_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:34:53) elapsed: 7.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/cutadapt/H9_PRO-seq_70.hist`

> `Degradation_ratio`	1.0311	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed_dups.fastq` (423122)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 4.435GB.  
  PID: 423122;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/processed_R2.flag` (423359)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 423359;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/repaired.flag`  

> `fastq_pair -t 72907421 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed.fastq` (423360)
<pre>
Left paired: 19677250		Right paired: 19677250
Left single: 156431		Right single: 1925340
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:20. Running peak memory: 4.435GB.  
  PID: 423360;	Command: fastq_pair;	Return code: 0;	Memory used: 3.441GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq` (423440)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 423440;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed.fastq` (423441)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 423441;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/repaired.flag` (423442)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 423442;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/dups_repaired.flag`  

> `fastq_pair -t 72907421 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed_dups.fastq` (423443)
<pre>
Left paired: 18110619		Right paired: 18110619
Left single: 116445		Right single: 3491971
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:20. Running peak memory: 4.435GB.  
  PID: 423443;	Command: fastq_pair;	Return code: 0;	Memory used: 3.305GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_trimmed.fastq` (423547)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 423547;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed_dups.fastq` (423549)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 423549;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/dups_repaired.flag` (423550)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 423550;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:37:58) elapsed: 185.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:37:58) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/human_rDNA_bt2` (423551)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 423551;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_R2.fq` (423552)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_70 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
19677250 reads; of these:
  19677250 (100.00%) were unpaired; of these:
    17637617 (89.63%) aligned 0 times
    2039633 (10.37%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.37% overall alignment rate

> `Aligned_reads_human_rDNA`	4079266.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.28	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:39:58) elapsed: 120.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/human_rDNA_dups_bt2` (423686)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 423686;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_dups_R2.fq` (423687)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_70 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/fastq/H9_PRO-seq_70_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2039633 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:41:52) elapsed: 114.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_70 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/tmpl0o9htnc -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp.bam` (424062,424063,424064)
<pre>
1636542 reads skipped
0 reads lost
17637617 reads; of these:
  17637617 (100.00%) were paired; of these:
    7312132 (41.46%) aligned concordantly 0 times
    8468481 (48.01%) aligned concordantly exactly 1 time
    1857004 (10.53%) aligned concordantly >1 times
    ----
    7312132 pairs aligned concordantly 0 times; of these:
      1626933 (22.25%) aligned discordantly 1 time
    ----
    5685199 pairs aligned 0 times concordantly or discordantly; of these:
      11370398 mates make up the pairs; of these:
        4666787 (41.04%) aligned 0 times
        2450638 (21.55%) aligned exactly 1 time
        4252973 (37.40%) aligned >1 times
86.77% overall alignment rate
[bam_sort_core] merging from 10 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:30:01. Running peak memory: 4.435GB.  
  PID: 424062;	Command: bowtie2;	Return code: 0;	Memory used: 3.77GB  
  PID: 424063;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 424064;	Command: samtools;	Return code: 0;	Memory used: 0.897GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort.bam` (427025)
<pre>
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 4.435GB.  
  PID: 427025;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	30608447	PEPPRO	_RES_

> `QC_filtered_reads`	17889191	PEPPRO	_RES_

> `Aligned_reads`	12719256.0	PEPPRO	_RES_

> `Alignment_rate`	32.06	PEPPRO	_RES_

> `Total_efficiency`	15.7	PEPPRO	_RES_

> `Read_depth`	2.64	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_70 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/tmpl0o9htnc -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp_dups.bam` (428504,428505,428506)
<pre>
16474077 reads; of these:
  16474077 (100.00%) were paired; of these:
    6685077 (40.58%) aligned concordantly 0 times
    8047929 (48.85%) aligned concordantly exactly 1 time
    1741071 (10.57%) aligned concordantly >1 times
    ----
    6685077 pairs aligned concordantly 0 times; of these:
      1550215 (23.19%) aligned discordantly 1 time
    ----
    5134862 pairs aligned 0 times concordantly or discordantly; of these:
      10269724 mates make up the pairs; of these:
        4278230 (41.66%) aligned 0 times
        2328888 (22.68%) aligned exactly 1 time
        3662606 (35.66%) aligned >1 times
87.02% overall alignment rate
[bam_sort_core] merging from 9 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:28:37. Running peak memory: 4.435GB.  
  PID: 428505;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 428504;	Command: bowtie2;	Return code: 0;	Memory used: 3.756GB  
  PID: 428506;	Command: samtools;	Return code: 0;	Memory used: 0.897GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort_dups.bam` (433052)
<pre>
</pre>
Command completed. Elapsed time: 0:00:57. Running peak memory: 4.435GB.  
  PID: 433052;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 08:52:27) elapsed: 4234.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_R2.fq` (433659)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 4.435GB.  
  PID: 433659;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/prealignments/H9_PRO-seq_70_human_rDNA_unmap_R1.fq` (433687)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 4.435GB.  
  PID: 433687;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp.bam` (433735)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 4.435GB.  
  PID: 433735;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	589030	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort.bam` (433805)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 4.435GB.  
  PID: 433805;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/chr_sizes.bed` (433821,433822,433823,433824)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 433823;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 433821;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 433824;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 433822;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_noMT.bam` (433826)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 4.435GB.  
  PID: 433826;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort.bam` (433870)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 433870;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort.bam` (433872)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 4.435GB.  
  PID: 433872;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Split BAM file (06-15 08:54:26) elapsed: 119.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam` (433890,433891)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:44. Running peak memory: 4.435GB.  
  PID: 433890;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 433891;	Command: samtools;	Return code: 0;	Memory used: 3.638GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE2.bam` (434298,434299)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:30. Running peak memory: 4.435GB.  
  PID: 434298;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 434299;	Command: samtools;	Return code: 0;	Memory used: 2.919GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp_dups.bam` (434564)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 4.435GB.  
  PID: 434564;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_dups_PE1.bam` (434635,434636)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 4.435GB.  
  PID: 434635;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 434636;	Command: samtools;	Return code: 0;	Memory used: 3.523GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_dups_PE2.bam` (435108,435110)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:28. Running peak memory: 4.435GB.  
  PID: 435108;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 435110;	Command: samtools;	Return code: 0;	Memory used: 2.827GB


### Calculate library complexity (06-15 09:01:55) elapsed: 449.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_dups_PE1.bam` (435289)
<pre>
BAM_INPUT
TOTAL READS     = 13202733
COUNTS_SUM      = 13202733
DISTINCT READS  = 1.0987e+07
DISTINCT COUNTS = 260
MAX COUNT       = 20256
COUNTS OF 1     = 9.91483e+06
OBSERVED COUNTS (20257)
1	9914831
2	738808
3	167534
4	65093
5	32981
6	18819
7	12049
8	7932
9	5463
10	4029
11	3004
12	2254
13	1773
14	1476
15	1227
16	1023
17	847
18	730
19	606
20	504
21	458
22	442
23	322
24	341
25	283
26	258
27	220
28	230
29	188
30	170
31	142
32	147
33	133
34	131
35	126
36	108
37	109
38	110
39	111
40	96
41	78
42	71
43	65
44	53
45	54
46	75
47	50
48	66
49	56
50	46
51	35
52	39
53	41
54	32
55	30
56	31
57	33
58	36
59	41
60	32
61	38
62	26
63	28
64	22
65	26
66	21
67	15
68	21
69	17
70	19
71	21
72	18
73	17
74	13
75	21
76	11
77	20
78	11
79	11
80	10
81	9
82	13
83	10
84	10
85	10
86	13
87	13
88	9
89	8
90	12
91	10
92	11
93	6
94	3
95	7
96	7
97	8
98	7
99	7
100	10
101	9
102	6
103	6
104	6
105	5
106	3
107	5
108	8
109	5
110	7
111	3
112	7
113	4
114	4
115	2
116	3
117	5
118	6
119	5
120	3
121	5
122	4
123	8
124	1
125	5
126	3
127	2
128	2
129	5
130	2
131	1
132	3
133	2
135	6
136	2
137	1
138	3
139	2
140	1
142	3
143	7
144	2
145	2
146	2
147	3
148	6
149	1
150	6
151	2
153	5
154	2
155	3
156	2
157	6
158	1
160	1
161	3
162	2
163	2
164	2
165	1
166	4
167	1
171	1
172	1
173	1
176	2
178	1
180	1
181	2
182	2
184	1
185	3
188	1
189	2
190	2
192	1
193	1
194	2
196	1
197	1
201	2
202	2
203	1
204	1
210	1
212	2
214	1
215	1
223	1
225	2
226	2
227	1
228	2
231	1
232	2
234	2
236	1
237	1
238	1
241	1
243	1
244	2
245	1
246	1
247	1
254	1
255	1
261	1
262	3
266	1
267	2
271	1
280	1
285	1
295	1
303	1
304	1
312	1
316	1
330	1
335	1
339	1
354	1
360	1
361	1
371	1
374	2
375	1
395	1
400	1
417	1
493	1
498	1
499	1
503	2
505	1
507	1
633	1
641	1
652	1
680	1
730	1
828	1
928	1
944	1
1153	1
1359	1
1363	1
1386	1
1655	1
1734	1
1854	1
2497	1
3223	1
6517	1
7685	1
9775	1
18952	1
20256	1

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
Command completed. Elapsed time: 0:01:12. Running peak memory: 4.435GB.  
  PID: 435289;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_dups_PE1.bam` (435463)
<pre>
BAM_INPUT
TOTAL READS     = 13202733
DISTINCT READS  = 1.0987e+07
DISTINCT COUNTS = 260
MAX COUNT       = 20256
COUNTS OF 1     = 9.91483e+06
MAX TERMS       = 100
OBSERVED COUNTS (20257)
1	9914831
2	738808
3	167534
4	65093
5	32981
6	18819
7	12049
8	7932
9	5463
10	4029
11	3004
12	2254
13	1773
14	1476
15	1227
16	1023
17	847
18	730
19	606
20	504
21	458
22	442
23	322
24	341
25	283
26	258
27	220
28	230
29	188
30	170
31	142
32	147
33	133
34	131
35	126
36	108
37	109
38	110
39	111
40	96
41	78
42	71
43	65
44	53
45	54
46	75
47	50
48	66
49	56
50	46
51	35
52	39
53	41
54	32
55	30
56	31
57	33
58	36
59	41
60	32
61	38
62	26
63	28
64	22
65	26
66	21
67	15
68	21
69	17
70	19
71	21
72	18
73	17
74	13
75	21
76	11
77	20
78	11
79	11
80	10
81	9
82	13
83	10
84	10
85	10
86	13
87	13
88	9
89	8
90	12
91	10
92	11
93	6
94	3
95	7
96	7
97	8
98	7
99	7
100	10
101	9
102	6
103	6
104	6
105	5
106	3
107	5
108	8
109	5
110	7
111	3
112	7
113	4
114	4
115	2
116	3
117	5
118	6
119	5
120	3
121	5
122	4
123	8
124	1
125	5
126	3
127	2
128	2
129	5
130	2
131	1
132	3
133	2
135	6
136	2
137	1
138	3
139	2
140	1
142	3
143	7
144	2
145	2
146	2
147	3
148	6
149	1
150	6
151	2
153	5
154	2
155	3
156	2
157	6
158	1
160	1
161	3
162	2
163	2
164	2
165	1
166	4
167	1
171	1
172	1
173	1
176	2
178	1
180	1
181	2
182	2
184	1
185	3
188	1
189	2
190	2
192	1
193	1
194	2
196	1
197	1
201	2
202	2
203	1
204	1
210	1
212	2
214	1
215	1
223	1
225	2
226	2
227	1
228	2
231	1
232	2
234	2
236	1
237	1
238	1
241	1
243	1
244	2
245	1
246	1
247	1
254	1
255	1
261	1
262	3
266	1
267	2
271	1
280	1
285	1
295	1
303	1
304	1
312	1
316	1
330	1
335	1
339	1
354	1
360	1
361	1
371	1
374	2
375	1
395	1
400	1
417	1
493	1
498	1
499	1
503	2
505	1
507	1
633	1
641	1
652	1
680	1
730	1
828	1
928	1
944	1
1153	1
1359	1
1363	1
1386	1
1655	1
1734	1
1854	1
2497	1
3223	1
6517	1
7685	1
9775	1
18952	1
20256	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
....................................................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:21. Running peak memory: 4.435GB.  
  PID: 435463;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_counts.txt` (435574)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 4.435GB.  
  PID: 435574;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_plot` (435613)
<pre>
Processing H9_PRO-seq_70
INFO: Found real counts for H9_PRO-seq_70 - Total (M): 13.641283 Unique (M): 13.202733

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.435GB.  
  PID: 435613;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Library complexity`	QC_hg38/H9_PRO-seq_70_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_70_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8537	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 09:04:53) elapsed: 178.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam` (435650)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.435GB.  
  PID: 435650;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_bamQC.tsv` (435859)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/tmp_H9_PRO-seq_70_PE1_th0jhiqu'
Processing with 12 cores...
Discarding 101 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270755v1']
Keeping 94 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 4.435GB.  
  PID: 435859;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.388GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	6820641.5	PEPPRO	_RES_

> `PBC2`	6820641.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_unmap.bam` (435909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.435GB.  
  PID: 435909;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_temp.bam`

> `Unmapped_reads`	4666787	PEPPRO	_RES_

### Split BAM by strand (06-15 09:05:26) elapsed: 33.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam` (435959)
<pre>
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 4.435GB.  
  PID: 435959;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam` (436486)
<pre>
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 4.435GB.  
  PID: 436486;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 09:07:01) elapsed: 96.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (436570)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 436570;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_plus_TssEnrichment.txt` (436583)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.435GB.  
  PID: 436583;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.981GB


> `TSS_coding_score`	33.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_minus_TssEnrichment.txt` (436634)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 4.435GB.  
  PID: 436634;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.103GB


> `TSS_non-coding_score`	12.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_minus_TssEnrichment.txt` (436668)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.435GB.  
  PID: 436668;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_70_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_70_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt` (436690,436691,436692,436693)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 436690;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 436692;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 436691;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 436693;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_keep.txt` (436695)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 436695;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 09:07:17) elapsed: 16.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_ensembl_tss.bed` (436697,436698)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 4.435GB.  
  PID: 436697;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 436698;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_ensembl_gene_body.bed` (436702,436703)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 436702;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 436703;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_TSS_density.bed` (436705,436706,436707,436708)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 4.435GB.  
  PID: 436705;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 436707;	Command: sort;	Return code: 0;	Memory used: 0.01GB  
  PID: 436706;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 436708;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_gene_body_density.bed` (436725,436726,436727)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 4.435GB.  
  PID: 436726;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 436725;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB  
  PID: 436727;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/tmpc1j1oqlj` (436757,436758,436759)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 436757;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 436759;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 436758;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/tmpc1j1oqlj | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0138162) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/tmpc1j1oqlj > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_pause_index.bed` (436765)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 436765;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	33.69	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_pause_index.bed` (436770)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.435GB.  
  PID: 436770;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/H9_PRO-seq_70_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_70_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_pause_index.bed` (436791)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 436791;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 09:08:06) elapsed: 49.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam`
12719256.0 4935879

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam`
12719256.0 4602962

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_gene_sort.bed` (436828,436829)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.435GB.  
  PID: 436828;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 436829;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_gene_coverage.bed` (436831)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 4.435GB.  
  PID: 436831;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/hg38_annotations.bed.gz` (436852)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 436852;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/hg38_annotations.bed` (436853)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 436853;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 09:08:53) elapsed: 47.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/raw/hg38_annotations.bed` (436862)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.435GB.  
  PID: 436862;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Enhancer_sort.bed` (436864,436865,436866,436867)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.435GB.  
  PID: 436864;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 436865;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 436867;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 436866;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Enhancer_plus_coverage.bed` (436869)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.435GB.  
  PID: 436869;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Enhancer_minus_coverage.bed` (436879)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.435GB.  
  PID: 436879;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_sort.bed` (436890,436891,436892,436893)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.435GB.  
  PID: 436890;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 436891;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 436893;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 436892;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_plus_coverage.bed` (436896)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 4.435GB.  
  PID: 436896;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_minus_coverage.bed` (436968)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.435GB.  
  PID: 436968;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_Flanking_Region"` (436977)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 436977;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed` (436978,436979,436980,436981)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.435GB.  
  PID: 436978;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 436980;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 436979;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 436981;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_Flanking_Region_plus_coverage.bed` (436984)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.435GB.  
  PID: 436984;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_Flanking_Region_minus_coverage.bed` (436993)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.435GB.  
  PID: 436993;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/5_UTR"` (437013)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 437013;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/5_UTR_sort.bed` (437014,437015,437016,437017)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.435GB.  
  PID: 437014;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 437015;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 437017;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 437016;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_5_UTR_plus_coverage.bed` (437019)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.435GB.  
  PID: 437019;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_5_UTR_minus_coverage.bed` (437290)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.435GB.  
  PID: 437290;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/3_UTR"` (437313)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 437313;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/3_UTR_sort.bed` (437314,437315,437316,437317)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.435GB.  
  PID: 437314;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 437315;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 437317;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB  
  PID: 437316;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_3_UTR_plus_coverage.bed` (437319)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.435GB.  
  PID: 437319;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_3_UTR_minus_coverage.bed` (437346)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 4.435GB.  
  PID: 437346;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Exon_sort.bed` (437356,437357,437358,437359)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 4.435GB.  
  PID: 437356;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 437357;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 437359;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 437358;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Exon_plus_coverage.bed` (437364)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 4.435GB.  
  PID: 437364;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Exon_minus_coverage.bed` (437383)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 4.435GB.  
  PID: 437383;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Intron_sort.bed` (437411,437412,437413,437414)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.435GB.  
  PID: 437411;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 437413;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 437412;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 437414;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Intron_plus_coverage.bed` (437417)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.435GB.  
  PID: 437417;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Intron_minus_coverage.bed` (437437)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.435GB.  
  PID: 437437;	Command: bedtools;	Return code: 0;	Memory used: 0.028GB


### Plot cFRiF/FRiF (06-15 09:11:17) elapsed: 144.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_70 -z 3099922541 -n 6962098 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Intron_plus_coverage.bed` (437472)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 4.435GB.  
  PID: 437472;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `cFRiF`	QC_hg38/H9_PRO-seq_70_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_70_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_70 -z 3099922541 -n 6962098 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_Intron_plus_coverage.bed` (437550)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 4.435GB.  
  PID: 437550;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `FRiF`	QC_hg38/H9_PRO-seq_70_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_70_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 09:12:19) elapsed: 62.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_exons_sort.bed` (438097,438098)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.435GB.  
  PID: 438098;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 438097;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_introns_sort.bed` (438104,438105,438106)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.435GB.  
  PID: 438104;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 438106;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 438105;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exons_coverage.bed` (438113)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 4.435GB.  
  PID: 438113;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_introns_coverage.bed` (438147)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 4.435GB.  
  PID: 438147;	Command: bedtools;	Return code: 0;	Memory used: 0.066GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/12.719256)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exons_rpkm.bed` (438184,438185,438186)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.435GB.  
  PID: 438184;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 438186;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 438185;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/12.719256)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_introns_rpkm.bed` (438189,438190,438191)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.435GB.  
  PID: 438189;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 438191;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 438190;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exon_intron_ratios.bed` (438193,438194,438195)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 438193;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 438195;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 438194;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exon_intron_ratios.bed --annotate` (438201)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.435GB.  
  PID: 438201;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_70_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_70_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/QC_hg38/H9_PRO-seq_70_exon_intron_ratios.bed` (438222)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.435GB.  
  PID: 438222;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 09:13:14) elapsed: 55.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam` (438230)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.435GB.  
  PID: 438230;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 12719256.0` (438236)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_70_plus_cuttrace_lwzndau9'
Processing with 4 cores...
stdin is empty of data
Discarding 115 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 80 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 80 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_plus_exact_body_0-mer.bw'
Merging 80 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:18. Running peak memory: 4.435GB.  
  PID: 438236;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.382GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam` (439527)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.435GB.  
  PID: 439527;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 12719256.0` (439532)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/aligned_hg38/H9_PRO-seq_70_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_70_minus_cuttrace_uvl4y49i'
Processing with 4 cores...
stdin is empty of data
Discarding 115 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270723v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 80 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 80 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_minus_exact_body_0-mer.bw'
Merging 80 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_70/signal_hg38/H9_PRO-seq_70_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:27. Running peak memory: 4.435GB.  
  PID: 439532;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.424GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:08:52
*  Total elapsed time (all runs):  3:43:48
*         Peak memory (this run):  4.4352 GB
*        Pipeline completed time: 2020-06-15 09:26:10
