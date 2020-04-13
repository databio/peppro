### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name H9_PRO-seq_3 --genome hg38 --input /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz --protocol PRO --umi-len 8 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj37-16c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/
*  Pipeline started at:   (02-27 09:26:03) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `16000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed/peppro/paper/dev4/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `H9_PRO-seq_3`
*              `scale`:  `False`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `paired`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `8`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz
Local input file: /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz

> `File_mb`	2340.86	PEPPRO	_RES_

> `Read_type`	paired	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-27 09:26:03) elapsed: 0.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz` (175991)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 175991;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz` (175992)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 175992;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq` (175993)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 0.003GB.  
  PID: 175993;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2.fastq` (393218)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 0.003GB.  
  PID: 393218;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	92628962	PEPPRO	_RES_

> `Fastq_reads`	92628962	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz']

### FASTQ processing:  (02-27 09:28:49) elapsed: 166.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt` (82667)
<pre>
</pre>
Command completed. Elapsed time: 0:01:02. Running peak memory: 0.284GB.  
  PID: 82667;	Command: cutadapt;	Return code: 0;	Memory used: 0.284GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (82751,82752)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 0.284GB.  
  PID: 82751;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 82752;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	27419170	PEPPRO	_RES_

> `Trim_loss_rate`	70.4	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (83025)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_3_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_3_R1_processed.fastq
Analysis complete for H9_PRO-seq_3_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 0.284GB.  
  PID: 83025;	Command: fastqc;	Return code: 0;	Memory used: 0.179GB

> `FastQC report r1`	fastqc/H9_PRO-seq_3_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_dedup.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq` (83116)
<pre>
[INFO][0m 1862033 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 1.955GB.  
  PID: 83116;	Command: seqkit;	Return code: 0;	Memory used: 1.955GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq` (83171,83172)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 1.955GB.  
  PID: 83171;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 83172;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	34725419.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	18106301.0	PEPPRO	_RES_

> `Duplicate_reads`	1862033.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	39.0943	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/processed_R1.flag` (83310)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.955GB.  
  PID: 83310;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R2_cutadapt.txt` (83313)
<pre>
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 1.955GB.  
  PID: 83313;	Command: cutadapt;	Return code: 0;	Memory used: 0.263GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (83401,83402)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 1.955GB.  
  PID: 83401;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 83402;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	54838340	PEPPRO	_RES_

> `Trim_loss_rate`	40.8	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (83625)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_3_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_3_R1_processed.fastq
Analysis complete for H9_PRO-seq_3_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 1.955GB.  
  PID: 83625;	Command: fastqc;	Return code: 0;	Memory used: 0.179GB

> `FastQC report r1`	fastqc/H9_PRO-seq_3_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (83713)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_3_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_3_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_3_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 1.955GB.  
  PID: 83713;	Command: fastqc;	Return code: 0;	Memory used: 0.17GB

> `FastQC report r2`	fastqc/H9_PRO-seq_3_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.histogram`  

> `fastq_pair -t 83366065 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq` (83804)
<pre>
Left paired: 28081513		Right paired: 28081513
Left single: 126667		Right single: 1382160
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:45. Running peak memory: 5.553GB.  
  PID: 83804;	Command: fastq_pair;	Return code: 0;	Memory used: 5.553GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq.paired.fq -o H9_PRO-seq_3 -d /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt` (84180)
<pre>
</pre>
Command completed. Elapsed time: 0:01:19. Running peak memory: 5.553GB.  
  PID: 84180;	Command: flash;	Return code: 0;	Memory used: 0.083GB


### Plot adapter insertion distribution (02-27 09:41:50) elapsed: 780.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt -u 8` (84412)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 5.553GB.  
  PID: 84412;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_3_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_3_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:41:55) elapsed: 6.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `Degradation_ratio`	0.9769	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq` (84439)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 5.553GB.  
  PID: 84439;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/processed_R2.flag` (84457)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.553GB.  
  PID: 84457;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/repaired.flag`  

> `fastq_pair -t 83366065 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (84458)
<pre>
Left paired: 27272197		Right paired: 27272197
Left single: 146973		Right single: 1431712
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:33. Running peak memory: 5.593GB.  
  PID: 84458;	Command: fastq_pair;	Return code: 0;	Memory used: 5.593GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (84616)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.593GB.  
  PID: 84616;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (84618)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.593GB.  
  PID: 84618;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/repaired.flag` (84619)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.593GB.  
  PID: 84619;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/dups_repaired.flag`  

> `fastq_pair -t 83366065 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq` (84620)
<pre>
Left paired: 25551175		Right paired: 25551175
Left single: 94450		Right single: 3152734
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:52. Running peak memory: 5.764GB.  
  PID: 84620;	Command: fastq_pair;	Return code: 0;	Memory used: 5.764GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq` (84980)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.764GB.  
  PID: 84980;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq` (84982)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 84982;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/dups_repaired.flag` (84985)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 84985;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-27 09:47:43) elapsed: 347.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:47:43) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2` (84986)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 84986;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq` (84987)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_3 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
27272197 reads; of these:
  27272197 (100.00%) were unpaired; of these:
    24569260 (90.09%) aligned 0 times
    2702937 (9.91%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.91% overall alignment rate

> `Aligned_reads_human_rDNA`	5405874.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.86	PEPPRO	_RES_

### Map to human_rDNA (02-27 09:50:57) elapsed: 194.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2` (85408)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 85408;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R2.fq` (85409)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_3 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2702937 reads skipped
0 reads lost

### Map to genome (02-27 09:54:03) elapsed: 187.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_3 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/tmpsl3pvknn -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam` (85608,85609,85610)
<pre>
2244415 reads skipped
0 reads lost
24569260 reads; of these:
  24569260 (100.00%) were paired; of these:
    9705626 (39.50%) aligned concordantly 0 times
    12263897 (49.92%) aligned concordantly exactly 1 time
    2599737 (10.58%) aligned concordantly >1 times
    ----
    9705626 pairs aligned concordantly 0 times; of these:
      2262313 (23.31%) aligned discordantly 1 time
    ----
    7443313 pairs aligned 0 times concordantly or discordantly; of these:
      14886626 mates make up the pairs; of these:
        5843431 (39.25%) aligned 0 times
        3306807 (22.21%) aligned exactly 1 time
        5736388 (38.53%) aligned >1 times
88.11% overall alignment rate
[bam_sort_core] merging from 13 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 1:02:37. Running peak memory: 5.764GB.  
  PID: 85609;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 85608;	Command: bowtie2;	Return code: 0;	Memory used: 3.661GB  
  PID: 85610;	Command: samtools;	Return code: 0;	Memory used: 0.874GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (92167)
<pre>
</pre>
Command completed. Elapsed time: 0:01:33. Running peak memory: 5.764GB.  
  PID: 92167;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	43295089	PEPPRO	_RES_

> `QC_filtered_reads`	25127415	PEPPRO	_RES_

> `Aligned_reads`	18167674.5	PEPPRO	_RES_

> `Alignment_rate`	33.13	PEPPRO	_RES_

> `Total_efficiency`	19.61	PEPPRO	_RES_

> `Read_depth`	2.98	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_3 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/tmpsl3pvknn -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam` (93786,93793,93795)
<pre>
23306760 reads; of these:
  23306760 (100.00%) were paired; of these:
    9052177 (38.84%) aligned concordantly 0 times
    11774634 (50.52%) aligned concordantly exactly 1 time
    2479949 (10.64%) aligned concordantly >1 times
    ----
    9052177 pairs aligned concordantly 0 times; of these:
      2177215 (24.05%) aligned discordantly 1 time
    ----
    6874962 pairs aligned 0 times concordantly or discordantly; of these:
      13749924 mates make up the pairs; of these:
        5423690 (39.45%) aligned 0 times
        3177751 (23.11%) aligned exactly 1 time
        5148483 (37.44%) aligned >1 times
88.36% overall alignment rate
[bam_sort_core] merging from 12 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:57:17. Running peak memory: 5.764GB.  
  PID: 93793;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 93786;	Command: bowtie2;	Return code: 0;	Memory used: 3.661GB  
  PID: 93795;	Command: samtools;	Return code: 0;	Memory used: 0.891GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam` (99610)
<pre>
</pre>
Command completed. Elapsed time: 0:01:36. Running peak memory: 5.764GB.  
  PID: 99610;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (02-27 12:10:39) elapsed: 8196.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq` (99942)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 5.764GB.  
  PID: 99942;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq` (99982)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 5.764GB.  
  PID: 99982;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam` (100035)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 5.764GB.  
  PID: 100035;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	711009	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (100071)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 5.764GB.  
  PID: 100071;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_noMT.bam` (100105,100106,100107,100108)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 5.764GB.  
  PID: 100105;	Command: samtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 100107;	Command: grep;	Return code: 0;	Memory used: 0.001GB  
  PID: 100106;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 100108;	Command: xargs;	Return code: 0;	Memory used: 0.059GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (100153)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 100153;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (100154)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 5.764GB.  
  PID: 100154;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


### Split BAM file (02-27 12:14:00) elapsed: 200.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam` (100188,100189)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:45. Running peak memory: 5.764GB.  
  PID: 100188;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 100189;	Command: samtools;	Return code: 0;	Memory used: 4.687GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE2.bam` (100575,100576)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:20. Running peak memory: 5.764GB.  
  PID: 100576;	Command: samtools;	Return code: 0;	Memory used: 3.727GB  
  PID: 100575;	Command: samtools;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam` (100795)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 5.764GB.  
  PID: 100795;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam` (101042,101043)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:52. Running peak memory: 5.764GB.  
  PID: 101042;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 101043;	Command: samtools;	Return code: 0;	Memory used: 4.574GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE2.bam` (101247,101248)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:24. Running peak memory: 5.764GB.  
  PID: 101247;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 101248;	Command: samtools;	Return code: 0;	Memory used: 3.715GB


### Calculate library complexity (02-27 12:25:40) elapsed: 700.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_out.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam` (101612)
<pre>
BAM_INPUT
TOTAL READS     = 19005606
COUNTS_SUM      = 19005606
DISTINCT READS  = 1.58447e+07
DISTINCT COUNTS = 285
MAX COUNT       = 19581
COUNTS OF 1     = 1.42036e+07
OBSERVED COUNTS (19582)
1	14203580
2	1156600
3	254232
4	94789
5	45836
6	26066
7	16009
8	10530
9	7388
10	5314
11	3937
12	3035
13	2237
14	1796
15	1517
16	1252
17	1031
18	892
19	767
20	636
21	535
22	462
23	411
24	386
25	367
26	297
27	272
28	227
29	236
30	215
31	209
32	162
33	160
34	160
35	150
36	126
37	130
38	108
39	116
40	107
41	94
42	95
43	81
44	70
45	69
46	67
47	68
48	59
49	70
50	59
51	60
52	65
53	42
54	43
55	50
56	44
57	45
58	37
59	39
60	31
61	47
62	32
63	28
64	33
65	38
66	28
67	31
68	33
69	29
70	31
71	21
72	20
73	18
74	21
75	24
76	18
77	28
78	20
79	17
80	17
81	14
82	17
83	14
84	13
85	13
86	17
87	11
88	17
89	14
90	8
91	13
92	10
93	15
94	14
95	11
96	5
97	13
98	10
99	10
100	12
101	3
102	14
103	11
104	8
105	15
106	9
107	9
108	9
109	6
110	9
111	5
112	4
113	4
114	6
115	12
116	1
117	7
118	6
119	9
120	4
121	2
122	3
123	4
124	5
125	4
126	3
127	4
128	7
129	7
130	2
131	5
132	3
133	5
134	8
135	3
136	5
137	4
138	7
139	2
140	2
141	3
142	4
143	3
144	6
145	3
146	6
147	2
148	1
149	5
150	2
151	3
152	4
153	5
154	2
155	5
156	3
157	4
158	5
159	2
160	4
161	2
162	1
163	3
164	1
166	3
167	1
168	3
169	3
170	3
172	3
173	6
174	1
175	2
176	2
177	1
178	2
179	1
181	1
183	1
184	3
185	3
186	1
187	1
188	2
189	1
190	1
191	2
193	3
194	1
195	2
196	2
197	3
199	3
200	2
201	2
202	2
203	4
204	1
205	1
207	4
208	1
210	1
211	1
212	1
213	2
214	2
216	1
217	2
223	3
224	1
225	1
226	1
229	2
231	2
233	1
234	1
235	1
238	1
239	1
242	2
246	1
247	1
248	3
250	1
254	1
256	2
259	2
260	1
264	1
271	1
274	1
276	1
284	2
294	1
299	1
304	1
305	1
318	1
334	1
335	1
348	1
351	1
353	1
357	1
366	1
390	2
391	1
393	1
394	1
397	1
422	1
432	1
443	1
450	1
467	1
469	1
500	1
508	1
513	1
516	1
526	1
537	1
543	1
560	1
577	1
582	1
611	1
619	1
653	1
697	1
887	1
939	1
964	1
1000	1
1280	1
1283	1
1718	1
1790	1
1929	1
2740	1
6648	1
7571	1
9539	1
16742	1
19581	1

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
sample size: 17000000
sample size: 18000000
sample size: 19000000
</pre>
Command completed. Elapsed time: 0:01:46. Running peak memory: 5.764GB.  
  PID: 101612;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam` (101715)
<pre>
BAM_INPUT
TOTAL READS     = 19005606
DISTINCT READS  = 1.58447e+07
DISTINCT COUNTS = 285
MAX COUNT       = 19581
COUNTS OF 1     = 1.42036e+07
MAX TERMS       = 100
OBSERVED COUNTS (19582)
1	14203580
2	1156600
3	254232
4	94789
5	45836
6	26066
7	16009
8	10530
9	7388
10	5314
11	3937
12	3035
13	2237
14	1796
15	1517
16	1252
17	1031
18	892
19	767
20	636
21	535
22	462
23	411
24	386
25	367
26	297
27	272
28	227
29	236
30	215
31	209
32	162
33	160
34	160
35	150
36	126
37	130
38	108
39	116
40	107
41	94
42	95
43	81
44	70
45	69
46	67
47	68
48	59
49	70
50	59
51	60
52	65
53	42
54	43
55	50
56	44
57	45
58	37
59	39
60	31
61	47
62	32
63	28
64	33
65	38
66	28
67	31
68	33
69	29
70	31
71	21
72	20
73	18
74	21
75	24
76	18
77	28
78	20
79	17
80	17
81	14
82	17
83	14
84	13
85	13
86	17
87	11
88	17
89	14
90	8
91	13
92	10
93	15
94	14
95	11
96	5
97	13
98	10
99	10
100	12
101	3
102	14
103	11
104	8
105	15
106	9
107	9
108	9
109	6
110	9
111	5
112	4
113	4
114	6
115	12
116	1
117	7
118	6
119	9
120	4
121	2
122	3
123	4
124	5
125	4
126	3
127	4
128	7
129	7
130	2
131	5
132	3
133	5
134	8
135	3
136	5
137	4
138	7
139	2
140	2
141	3
142	4
143	3
144	6
145	3
146	6
147	2
148	1
149	5
150	2
151	3
152	4
153	5
154	2
155	5
156	3
157	4
158	5
159	2
160	4
161	2
162	1
163	3
164	1
166	3
167	1
168	3
169	3
170	3
172	3
173	6
174	1
175	2
176	2
177	1
178	2
179	1
181	1
183	1
184	3
185	3
186	1
187	1
188	2
189	1
190	1
191	2
193	3
194	1
195	2
196	2
197	3
199	3
200	2
201	2
202	2
203	4
204	1
205	1
207	4
208	1
210	1
211	1
212	1
213	2
214	2
216	1
217	2
223	3
224	1
225	1
226	1
229	2
231	2
233	1
234	1
235	1
238	1
239	1
242	2
246	1
247	1
248	3
250	1
254	1
256	2
259	2
260	1
264	1
271	1
274	1
276	1
284	2
294	1
299	1
304	1
305	1
318	1
334	1
335	1
348	1
351	1
353	1
357	1
366	1
390	2
391	1
393	1
394	1
397	1
422	1
432	1
443	1
450	1
467	1
469	1
500	1
508	1
513	1
516	1
526	1
537	1
543	1
560	1
577	1
582	1
611	1
619	1
653	1
697	1
887	1
939	1
964	1
1000	1
1280	1
1283	1
1718	1
1790	1
1929	1
2740	1
6648	1
7571	1
9539	1
16742	1
19581	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
....................................................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:54. Running peak memory: 5.764GB.  
  PID: 101715;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_counts.txt` (101830)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 5.764GB.  
  PID: 101830;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_plot` (101871)
<pre>
Processing H9_PRO-seq_3
INFO: Found real counts for H9_PRO-seq_3 - Total (M): 19.481041 Unique (M): 19.005606

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.764GB.  
  PID: 101871;	Command: Rscript;	Return code: 0;	Memory used: 0.237GB

> `Library complexity`	QC_hg38/H9_PRO-seq_3_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_3_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8825	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 12:29:54) elapsed: 254.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam` (101890)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 5.764GB.  
  PID: 101890;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv` (102130)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/tmp_H9_PRO-seq_3_PE1_dkf5d0_7'
Processing with 8 cores...
Discarding 101 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 94 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 5.764GB.  
  PID: 102130;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.507GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	9740520.5	PEPPRO	_RES_

> `PBC2`	9740520.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_unmap.bam` (102172)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 5.764GB.  
  PID: 102172;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam`

> `Unmapped_reads`	5843431	PEPPRO	_RES_

### Split BAM by strand (02-27 12:30:48) elapsed: 54.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam` (102220)
<pre>
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 5.764GB.  
  PID: 102220;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam` (102294)
<pre>
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 5.764GB.  
  PID: 102294;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 12:33:10) elapsed: 142.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (102368)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 102368;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt` (102370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 5.764GB.  
  PID: 102370;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.489GB


> `TSS_coding_score`	28.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt` (102427)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 5.764GB.  
  PID: 102427;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.486GB


> `TSS_non-coding_score`	10.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt` (102452)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 5.764GB.  
  PID: 102452;	Command: Rscript;	Return code: 0;	Memory used: 0.335GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_3_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_3_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt` (102470,102471,102472,102473)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 102470;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 102472;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 102471;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 102473;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt` (102475)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 102475;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (02-27 12:33:30) elapsed: 20.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed` (102477,102478)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 5.764GB.  
  PID: 102477;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 102478;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed` (102482,102483)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 102482;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 102483;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSS_density.bed` (102485,102486,102487,102488)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 5.764GB.  
  PID: 102486;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 102488;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 102485;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 102487;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_gene_body_density.bed` (102523,102524,102525)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 5.764GB.  
  PID: 102525;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 102523;	Command: bedtools;	Return code: 0;	Memory used: 0.062GB  
  PID: 102524;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed` (102554,102555,102556)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 102554;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 102556;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 102555;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	33.89	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed` (102561)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 5.764GB.  
  PID: 102561;	Command: Rscript;	Return code: 0;	Memory used: 0.366GB

> `Pause index`	QC_hg38/H9_PRO-seq_3_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_3_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed` (102584)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 102584;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 12:34:36) elapsed: 67.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam`
18167674.5 7066591

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam`
18167674.5 6601482

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_gene_sort.bed` (102832,102833)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.764GB.  
  PID: 102832;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 102833;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_gene_coverage.bed` (102836)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 5.764GB.  
  PID: 102836;	Command: bedtools;	Return code: 0;	Memory used: 0.067GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed.gz` (102874)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 102874;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed` (102875)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 102875;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 12:35:46) elapsed: 70.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed` (102884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.764GB.  
  PID: 102884;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR"` (102887)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 102887;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed` (102888,102889,102890,102891)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.764GB.  
  PID: 102888;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 102889;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 102891;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 102890;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed` (102894)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.764GB.  
  PID: 102894;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_minus_coverage.bed` (102907)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.764GB.  
  PID: 102907;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR"` (102919)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 102919;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed` (102920,102921,102922,102923)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.764GB.  
  PID: 102920;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 102921;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 102923;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 102922;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed` (102926)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 5.764GB.  
  PID: 102926;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_minus_coverage.bed` (102939)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.764GB.  
  PID: 102939;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed` (102960,102961,102962,102963)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.764GB.  
  PID: 102960;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 102961;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 102963;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB  
  PID: 102962;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed` (102965)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.764GB.  
  PID: 102965;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_minus_coverage.bed` (102981)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.764GB.  
  PID: 102981;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed` (102993,102994,102995,102996)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 5.764GB.  
  PID: 102993;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 102994;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 102996;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 102995;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed` (103001)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 5.764GB.  
  PID: 103001;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_minus_coverage.bed` (103014)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 5.764GB.  
  PID: 103014;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed` (103035,103036,103037,103038)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.764GB.  
  PID: 103035;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 103037;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 103036;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 103038;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed` (103042)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 5.764GB.  
  PID: 103042;	Command: bedtools;	Return code: 0;	Memory used: 0.06GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_minus_coverage.bed` (103059)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 5.764GB.  
  PID: 103059;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed` (103075,103076,103077,103078)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 103075;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 103076;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 103078;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 103077;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed` (103080)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 5.764GB.  
  PID: 103080;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_minus_coverage.bed` (103093)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.764GB.  
  PID: 103093;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region"` (103113)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 103113;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed` (103114,103115,103116,103117)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.764GB.  
  PID: 103114;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 103116;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 103115;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 103117;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed` (103120)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.764GB.  
  PID: 103120;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_minus_coverage.bed` (103136)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.764GB.  
  PID: 103136;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


### Plot cFRiF/FRiF (02-27 12:39:09) elapsed: 203.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_3 -z 3099922541 -n 9939223 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed` (103160)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 5.764GB.  
  PID: 103160;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `cFRiF`	QC_hg38/H9_PRO-seq_3_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_3_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_3 -z 3099922541 -n 9939223 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed` (103208)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 5.764GB.  
  PID: 103208;	Command: Rscript;	Return code: 0;	Memory used: 0.465GB

> `FRiF`	QC_hg38/H9_PRO-seq_3_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_3_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 12:40:09) elapsed: 60.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_exons_sort.bed` (103453,103454)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.764GB.  
  PID: 103454;	Command: bedtools;	Return code: 0;	Memory used: 0.085GB  
  PID: 103453;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_introns_sort.bed` (103464,103465,103466)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 5.764GB.  
  PID: 103464;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 103466;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 103465;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_coverage.bed` (103473)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 5.764GB.  
  PID: 103473;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_coverage.bed` (103505)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 5.764GB.  
  PID: 103505;	Command: bedtools;	Return code: 0;	Memory used: 0.042GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.1676745)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_rpkm.bed` (103533,103535,103536)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.764GB.  
  PID: 103533;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 103536;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 103535;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.1676745)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_rpkm.bed` (103538,103539,103540)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.764GB.  
  PID: 103538;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 103540;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 103539;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed` (103543,103544,103545)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 103543;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 103545;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 103544;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.28	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed --annotate` (103551)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 5.764GB.  
  PID: 103551;	Command: Rscript;	Return code: 0;	Memory used: 0.121GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_3_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_3_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed` (103571)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.764GB.  
  PID: 103571;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-27 12:41:25) elapsed: 76.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam` (103580)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 5.764GB.  
  PID: 103580;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (103587)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_3_plus_cuttrace_6dkzkl3l'
Processing with 2 cores...
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:22. Running peak memory: 5.764GB.  
  PID: 103587;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.428GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam` (105396)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 5.764GB.  
  PID: 105396;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (105403)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_3_minus_cuttrace_5hh3rwrj'
Processing with 2 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_exact_body_0-mer.bw'
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:17. Running peak memory: 5.764GB.  
  PID: 105403;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.758GB

Starting cleanup: 77 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:36:20
*  Total elapsed time (all runs):  7:04:55
*         Peak memory (this run):  5.7636 GB
*        Pipeline completed time: 2020-02-27 13:02:22
