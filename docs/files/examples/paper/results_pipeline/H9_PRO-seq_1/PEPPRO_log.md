### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name H9_PRO-seq_1 --genome hg38 --input /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz --protocol PRO --umi-len 8 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj37-15c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/
*  Pipeline started at:   (02-27 09:25:33) elapsed: 1.0 _TIME_

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
*              `input`:  `['/project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_1`
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

Local input file: /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz
Local input file: /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz

> `File_mb`	2469.73	PEPPRO	_RES_

> `Read_type`	paired	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-27 09:25:33) elapsed: 0.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz` (161245)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 161245;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz` (161246)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 161246;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq` (161247)
<pre>
</pre>
Command completed. Elapsed time: 0:00:52. Running peak memory: 0.003GB.  
  PID: 161247;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2.fastq` (389922)
<pre>
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 0.003GB.  
  PID: 389922;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	97729860	PEPPRO	_RES_

> `Fastq_reads`	97729860	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz']

### FASTQ processing:  (02-27 09:28:04) elapsed: 151.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt` (112262)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 0.277GB.  
  PID: 112262;	Command: cutadapt;	Return code: 0;	Memory used: 0.277GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (112336,112337)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 0.277GB.  
  PID: 112336;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 112337;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	30221343	PEPPRO	_RES_

> `Trim_loss_rate`	69.08	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (112370)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_1_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_1_R1_processed.fastq
Analysis complete for H9_PRO-seq_1_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 0.277GB.  
  PID: 112370;	Command: fastqc;	Return code: 0;	Memory used: 0.18GB

> `FastQC report r1`	fastqc/H9_PRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_dedup.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq` (112681)
<pre>
[INFO][0m 3516802 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 1.015GB.  
  PID: 112681;	Command: seqkit;	Return code: 0;	Memory used: 1.015GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq` (112727,112728)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 1.015GB.  
  PID: 112727;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 112728;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	34822522.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	17998987.0	PEPPRO	_RES_

> `Duplicate_reads`	3516802.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	36.8342	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/processed_R1.flag` (112793)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.015GB.  
  PID: 112793;	Command: touch;	Return code: 0;	Memory used: 0.002GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R2_cutadapt.txt` (112796)
<pre>
</pre>
Command completed. Elapsed time: 0:01:02. Running peak memory: 1.015GB.  
  PID: 112796;	Command: cutadapt;	Return code: 0;	Memory used: 0.274GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (112905,112906)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 1.015GB.  
  PID: 112905;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 112906;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	60442686	PEPPRO	_RES_

> `Trim_loss_rate`	38.15	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (112941)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_1_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_1_R1_processed.fastq
Analysis complete for H9_PRO-seq_1_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:07. Running peak memory: 1.015GB.  
  PID: 112941;	Command: fastqc;	Return code: 0;	Memory used: 0.18GB

> `FastQC report r1`	fastqc/H9_PRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (113016)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_1_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_1_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_1_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:01:07. Running peak memory: 1.015GB.  
  PID: 113016;	Command: fastqc;	Return code: 0;	Memory used: 0.172GB

> `FastQC report r2`	fastqc/H9_PRO-seq_1_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.histogram`  

> `fastq_pair -t 87956874 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq` (113303)
<pre>
Left paired: 30692275		Right paired: 30692275
Left single: 173668		Right single: 1557501
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:31. Running peak memory: 6.358GB.  
  PID: 113303;	Command: fastq_pair;	Return code: 0;	Memory used: 6.358GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq.paired.fq -o H9_PRO-seq_1 -d /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt` (113436)
<pre>
</pre>
Command completed. Elapsed time: 0:01:21. Running peak memory: 6.358GB.  
  PID: 113436;	Command: flash;	Return code: 0;	Memory used: 0.122GB


### Plot adapter insertion distribution (02-27 09:39:50) elapsed: 705.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt -u 8` (113647)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.358GB.  
  PID: 113647;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:39:56) elapsed: 7.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `Degradation_ratio`	0.8896	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq` (113678)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 6.358GB.  
  PID: 113678;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/processed_R2.flag` (113917)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.358GB.  
  PID: 113917;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/repaired.flag`  

> `fastq_pair -t 87956874 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (113918)
<pre>
Left paired: 30003416		Right paired: 30003416
Left single: 217927		Right single: 1602455
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:12. Running peak memory: 6.382GB.  
  PID: 113918;	Command: fastq_pair;	Return code: 0;	Memory used: 6.382GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (114031)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.382GB.  
  PID: 114031;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (114033)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 114033;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/repaired.flag` (114034)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 114034;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/dups_repaired.flag`  

> `fastq_pair -t 87956874 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq` (114035)
<pre>
Left paired: 26674489		Right paired: 26674489
Left single: 126830		Right single: 4931382
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:01. Running peak memory: 6.382GB.  
  PID: 114035;	Command: fastq_pair;	Return code: 0;	Memory used: 6.017GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq` (114141)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.382GB.  
  PID: 114141;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq` (114143)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 114143;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/dups_repaired.flag` (114145)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 114145;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-27 09:44:32) elapsed: 275.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:44:32) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2` (114146)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 114146;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq` (114147)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
30003416 reads; of these:
  30003416 (100.00%) were unpaired; of these:
    27299076 (90.99%) aligned 0 times
    2704340 (9.01%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.01% overall alignment rate

> `Aligned_reads_human_rDNA`	5408680.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	8.95	PEPPRO	_RES_

### Map to human_rDNA (02-27 09:47:36) elapsed: 184.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2` (114526)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 114526;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R2.fq` (114527)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2704340 reads skipped
0 reads lost

### Map to genome (02-27 09:50:18) elapsed: 162.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_1 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/tmpxx1gzixa -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam` (114903,114904,114906)
<pre>
2176982 reads skipped
0 reads lost
27299076 reads; of these:
  27299076 (100.00%) were paired; of these:
    9566070 (35.04%) aligned concordantly 0 times
    14769056 (54.10%) aligned concordantly exactly 1 time
    2963950 (10.86%) aligned concordantly >1 times
    ----
    9566070 pairs aligned concordantly 0 times; of these:
      2129010 (22.26%) aligned discordantly 1 time
    ----
    7437060 pairs aligned 0 times concordantly or discordantly; of these:
      14874120 mates make up the pairs; of these:
        6258301 (42.08%) aligned 0 times
        3282618 (22.07%) aligned exactly 1 time
        5333201 (35.86%) aligned >1 times
88.54% overall alignment rate
[bam_sort_core] merging from 14 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 1:01:04. Running peak memory: 6.382GB.  
  PID: 114903;	Command: bowtie2;	Return code: 0;	Memory used: 3.698GB  
  PID: 114904;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 114906;	Command: samtools;	Return code: 0;	Memory used: 0.877GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (120795)
<pre>
</pre>
Command completed. Elapsed time: 0:01:29. Running peak memory: 6.382GB.  
  PID: 120795;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	48339851	PEPPRO	_RES_

> `QC_filtered_reads`	27585022	PEPPRO	_RES_

> `Aligned_reads`	20754829.0	PEPPRO	_RES_

> `Alignment_rate`	34.34	PEPPRO	_RES_

> `Total_efficiency`	21.24	PEPPRO	_RES_

> `Read_depth`	3.26	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_1 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/tmpxx1gzixa -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam` (122293,122299,122300)
<pre>
24497507 reads; of these:
  24497507 (100.00%) were paired; of these:
    8334941 (34.02%) aligned concordantly 0 times
    13482466 (55.04%) aligned concordantly exactly 1 time
    2680100 (10.94%) aligned concordantly >1 times
    ----
    8334941 pairs aligned concordantly 0 times; of these:
      1950931 (23.41%) aligned discordantly 1 time
    ----
    6384010 pairs aligned 0 times concordantly or discordantly; of these:
      12768020 mates make up the pairs; of these:
        5376104 (42.11%) aligned 0 times
        3000379 (23.50%) aligned exactly 1 time
        4391537 (34.39%) aligned >1 times
89.03% overall alignment rate
[bam_sort_core] merging from 13 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:54:06. Running peak memory: 6.382GB.  
  PID: 122299;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 122293;	Command: bowtie2;	Return code: 0;	Memory used: 3.664GB  
  PID: 122300;	Command: samtools;	Return code: 0;	Memory used: 0.873GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam` (127587)
<pre>
</pre>
Command completed. Elapsed time: 0:01:19. Running peak memory: 6.382GB.  
  PID: 127587;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (02-27 12:01:42) elapsed: 7884.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq` (127683)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 6.382GB.  
  PID: 127683;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq` (127724)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.382GB.  
  PID: 127724;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam` (127794)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 6.382GB.  
  PID: 127794;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	674956	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (127833)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.382GB.  
  PID: 127833;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_noMT.bam` (127857,127858,127859,127860)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 6.382GB.  
  PID: 127857;	Command: samtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 127859;	Command: grep;	Return code: 0;	Memory used: 0.001GB  
  PID: 127858;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 127860;	Command: xargs;	Return code: 0;	Memory used: 0.055GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (127908)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 127908;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (127909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 6.382GB.  
  PID: 127909;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


### Split BAM file (02-27 12:05:00) elapsed: 198.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam` (127937,128061)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:47. Running peak memory: 6.382GB.  
  PID: 127937;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 128061;	Command: samtools;	Return code: 0;	Memory used: 5.306GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE2.bam` (128319,128320)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:37. Running peak memory: 6.382GB.  
  PID: 128319;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 128320;	Command: samtools;	Return code: 0;	Memory used: 4.424GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam` (128760)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 6.382GB.  
  PID: 128760;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam` (128823,128824)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:38. Running peak memory: 6.382GB.  
  PID: 128823;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 128824;	Command: samtools;	Return code: 0;	Memory used: 4.914GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE2.bam` (128994,128995)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:19. Running peak memory: 6.382GB.  
  PID: 128995;	Command: samtools;	Return code: 0;	Memory used: 4.023GB  
  PID: 128994;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


### Calculate library complexity (02-27 12:16:35) elapsed: 695.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_out.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam` (129390)
<pre>
BAM_INPUT
TOTAL READS     = 20330230
COUNTS_SUM      = 20330230
DISTINCT READS  = 1.66542e+07
DISTINCT COUNTS = 262
MAX COUNT       = 16193
COUNTS OF 1     = 1.47277e+07
OBSERVED COUNTS (16194)
1	14727710
2	1349104
3	301309
4	112637
5	55354
6	31165
7	18956
8	12770
9	8912
10	6479
11	4906
12	3639
13	2910
14	2335
15	1909
16	1594
17	1295
18	1055
19	904
20	791
21	668
22	611
23	529
24	477
25	427
26	388
27	364
28	323
29	242
30	234
31	263
32	232
33	207
34	182
35	165
36	158
37	167
38	119
39	133
40	114
41	111
42	122
43	105
44	82
45	73
46	82
47	75
48	75
49	76
50	60
51	49
52	59
53	63
54	55
55	46
56	48
57	54
58	37
59	41
60	42
61	47
62	48
63	35
64	35
65	32
66	27
67	30
68	32
69	27
70	18
71	21
72	23
73	16
74	16
75	21
76	17
77	15
78	13
79	24
80	11
81	16
82	21
83	21
84	11
85	11
86	10
87	23
88	12
89	16
90	16
91	11
92	12
93	10
94	8
95	9
96	12
97	16
98	6
99	14
100	3
101	5
102	7
103	9
104	5
105	11
106	10
107	7
108	5
109	7
110	7
111	7
112	4
113	7
114	3
115	4
116	11
117	9
118	4
119	5
122	7
123	3
124	3
125	8
126	3
127	4
128	4
129	3
130	3
131	7
132	6
133	4
134	4
135	7
136	7
137	2
138	5
139	4
140	3
141	1
142	4
143	4
144	2
145	3
146	4
147	5
148	2
149	4
150	1
151	1
152	1
154	1
155	5
156	3
158	1
161	1
162	3
163	2
164	1
165	3
166	3
167	1
169	1
170	2
171	1
172	1
173	2
174	1
175	1
176	3
178	2
179	4
181	1
182	1
183	1
184	1
186	1
187	2
188	1
190	1
192	2
193	1
196	1
198	1
199	1
200	1
202	1
203	2
204	2
205	1
206	1
207	1
208	1
209	1
212	1
214	1
216	3
220	1
221	1
222	1
223	1
224	2
227	1
228	1
232	1
233	1
235	1
237	2
241	1
249	1
251	1
253	1
254	1
255	1
256	1
257	1
258	2
262	1
264	1
266	1
268	1
269	1
275	1
277	1
280	1
296	1
299	1
303	1
306	1
313	1
315	1
317	2
325	1
345	1
348	1
406	1
412	1
430	1
433	1
436	1
437	1
493	1
494	1
540	1
547	1
606	2
611	1
612	1
734	1
787	1
789	1
848	1
927	1
943	1
1376	1
1390	1
1414	1
2008	1
4768	1
5683	1
6998	1
14724	1
16193	1

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
sample size: 20000000
</pre>
Command completed. Elapsed time: 0:01:41. Running peak memory: 6.382GB.  
  PID: 129390;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam` (129478)
<pre>
BAM_INPUT
TOTAL READS     = 20330230
DISTINCT READS  = 1.66542e+07
DISTINCT COUNTS = 262
MAX COUNT       = 16193
COUNTS OF 1     = 1.47277e+07
MAX TERMS       = 100
OBSERVED COUNTS (16194)
1	14727710
2	1349104
3	301309
4	112637
5	55354
6	31165
7	18956
8	12770
9	8912
10	6479
11	4906
12	3639
13	2910
14	2335
15	1909
16	1594
17	1295
18	1055
19	904
20	791
21	668
22	611
23	529
24	477
25	427
26	388
27	364
28	323
29	242
30	234
31	263
32	232
33	207
34	182
35	165
36	158
37	167
38	119
39	133
40	114
41	111
42	122
43	105
44	82
45	73
46	82
47	75
48	75
49	76
50	60
51	49
52	59
53	63
54	55
55	46
56	48
57	54
58	37
59	41
60	42
61	47
62	48
63	35
64	35
65	32
66	27
67	30
68	32
69	27
70	18
71	21
72	23
73	16
74	16
75	21
76	17
77	15
78	13
79	24
80	11
81	16
82	21
83	21
84	11
85	11
86	10
87	23
88	12
89	16
90	16
91	11
92	12
93	10
94	8
95	9
96	12
97	16
98	6
99	14
100	3
101	5
102	7
103	9
104	5
105	11
106	10
107	7
108	5
109	7
110	7
111	7
112	4
113	7
114	3
115	4
116	11
117	9
118	4
119	5
122	7
123	3
124	3
125	8
126	3
127	4
128	4
129	3
130	3
131	7
132	6
133	4
134	4
135	7
136	7
137	2
138	5
139	4
140	3
141	1
142	4
143	4
144	2
145	3
146	4
147	5
148	2
149	4
150	1
151	1
152	1
154	1
155	5
156	3
158	1
161	1
162	3
163	2
164	1
165	3
166	3
167	1
169	1
170	2
171	1
172	1
173	2
174	1
175	1
176	3
178	2
179	4
181	1
182	1
183	1
184	1
186	1
187	2
188	1
190	1
192	2
193	1
196	1
198	1
199	1
200	1
202	1
203	2
204	2
205	1
206	1
207	1
208	1
209	1
212	1
214	1
216	3
220	1
221	1
222	1
223	1
224	2
227	1
228	1
232	1
233	1
235	1
237	2
241	1
249	1
251	1
253	1
254	1
255	1
256	1
257	1
258	2
262	1
264	1
266	1
268	1
269	1
275	1
277	1
280	1
296	1
299	1
303	1
306	1
313	1
315	1
317	2
325	1
345	1
348	1
406	1
412	1
430	1
433	1
436	1
437	1
493	1
494	1
540	1
547	1
606	2
611	1
612	1
734	1
787	1
789	1
848	1
927	1
943	1
1376	1
1390	1
1414	1
2008	1
4768	1
5683	1
6998	1
14724	1
16193	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
....................................................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:49. Running peak memory: 6.382GB.  
  PID: 129478;	Command: preseq;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_counts.txt` (129800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.382GB.  
  PID: 129800;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_plot` (129832)
<pre>
Processing H9_PRO-seq_1
INFO: Found real counts for H9_PRO-seq_1 - Total (M): 21.966708 Unique (M): 20.33023

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 6.382GB.  
  PID: 129832;	Command: Rscript;	Return code: 0;	Memory used: 0.367GB

> `Library complexity`	QC_hg38/H9_PRO-seq_1_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_1_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8774	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 12:20:38) elapsed: 242.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam` (129850)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.382GB.  
  PID: 129850;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv` (129866)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/tmp_H9_PRO-seq_1_PE1_bn8507ak'
Processing with 8 cores...
Discarding 96 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270720v1_random', 'chr14_KI270724v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 99 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270539v1', 'chrUn_KI270587v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 6.382GB.  
  PID: 129866;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.916GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	10983354.0	PEPPRO	_RES_

> `PBC2`	10983354.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_unmap.bam` (129909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.382GB.  
  PID: 129909;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam`

> `Unmapped_reads`	6258301	PEPPRO	_RES_

### Split BAM by strand (02-27 12:21:34) elapsed: 57.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam` (129947)
<pre>
</pre>
Command completed. Elapsed time: 0:01:11. Running peak memory: 6.382GB.  
  PID: 129947;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam` (130010)
<pre>
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 6.382GB.  
  PID: 130010;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 12:23:54) elapsed: 139.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (130074)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130074;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt` (130075)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.382GB.  
  PID: 130075;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.607GB


> `TSS_coding_score`	26.7	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt` (130100)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.382GB.  
  PID: 130100;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.624GB


> `TSS_non-coding_score`	10.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt` (130126)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 6.382GB.  
  PID: 130126;	Command: Rscript;	Return code: 0;	Memory used: 0.343GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_1_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_1_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt` (130144,130145,130146,130147)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130144;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 130146;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 130145;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 130147;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt` (130149)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130149;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-27 12:24:13) elapsed: 20.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed` (130151,130152)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.382GB.  
  PID: 130151;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 130152;	Command: bedtools;	Return code: 0;	Memory used: 0.1GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed` (130156,130157)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130156;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 130157;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSS_density.bed` (130159,130160,130161,130162)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.382GB.  
  PID: 130160;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 130162;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 130159;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 130161;	Command: sort;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_gene_body_density.bed` (130186,130187,130188)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.382GB.  
  PID: 130187;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 130186;	Command: bedtools;	Return code: 0;	Memory used: 0.063GB  
  PID: 130188;	Command: sort;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed` (130415,130416,130417)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130415;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 130417;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 130416;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	32.63	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed` (130422)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.382GB.  
  PID: 130422;	Command: Rscript;	Return code: 0;	Memory used: 0.238GB

> `Pause index`	QC_hg38/H9_PRO-seq_1_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_1_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed` (130442)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130442;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 12:25:19) elapsed: 66.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam`
20754829.0 7959339

> `Plus_FRiP`	0.38	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam`
20754829.0 7448642

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_gene_sort.bed` (130506,130507)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.382GB.  
  PID: 130506;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 130507;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_gene_coverage.bed` (130510)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 6.382GB.  
  PID: 130510;	Command: bedtools;	Return code: 0;	Memory used: 0.069GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed.gz` (130537)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130537;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed` (130538)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.382GB.  
  PID: 130538;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 12:26:30) elapsed: 71.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed` (130547)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.382GB.  
  PID: 130547;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR"` (130549)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130549;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed` (130550,130551,130552,130553)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.382GB.  
  PID: 130550;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 130551;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 130553;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 130552;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed` (130557)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.382GB.  
  PID: 130557;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_minus_coverage.bed` (130570)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.382GB.  
  PID: 130570;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR"` (130585)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130585;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed` (130586,130587,130588,130589)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.382GB.  
  PID: 130586;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 130587;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 130589;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 130588;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed` (130591)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.382GB.  
  PID: 130591;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_minus_coverage.bed` (130604)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.382GB.  
  PID: 130604;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed` (130616,130617,130618,130619)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.382GB.  
  PID: 130616;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 130617;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 130619;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB  
  PID: 130618;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed` (130621)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.382GB.  
  PID: 130621;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_minus_coverage.bed` (130634)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.382GB.  
  PID: 130634;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed` (130647,130648,130649,130650)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.382GB.  
  PID: 130647;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 130648;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 130650;	Command: bedtools;	Return code: 0;	Memory used: 0.18GB  
  PID: 130649;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed` (130658)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.382GB.  
  PID: 130658;	Command: bedtools;	Return code: 0;	Memory used: 0.026GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_minus_coverage.bed` (130671)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.382GB.  
  PID: 130671;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed` (130684,130685,130686,130687)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.382GB.  
  PID: 130684;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 130686;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 130685;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 130687;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed` (130690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.382GB.  
  PID: 130690;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_minus_coverage.bed` (130705)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.382GB.  
  PID: 130705;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed` (130722,130723,130724,130725)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130722;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 130723;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 130725;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 130724;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed` (130727)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.382GB.  
  PID: 130727;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_minus_coverage.bed` (130741)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.382GB.  
  PID: 130741;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region"` (130754)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 130754;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed` (130755,130756,130757,130758)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.382GB.  
  PID: 130755;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 130757;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 130756;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 130758;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (130761)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.382GB.  
  PID: 130761;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed` (130774)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.382GB.  
  PID: 130774;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


### Plot cFRiF/FRiF (02-27 12:29:50) elapsed: 200.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_1 -z 3099922541 -n 11202246 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (130797)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 6.382GB.  
  PID: 130797;	Command: Rscript;	Return code: 0;	Memory used: 0.468GB

> `cFRiF`	QC_hg38/H9_PRO-seq_1_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_1_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_1 -z 3099922541 -n 11202246 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (131072)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 6.382GB.  
  PID: 131072;	Command: Rscript;	Return code: 0;	Memory used: 0.47GB

> `FRiF`	QC_hg38/H9_PRO-seq_1_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_1_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 12:30:44) elapsed: 54.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_exons_sort.bed` (131102,131103)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.382GB.  
  PID: 131103;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB  
  PID: 131102;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_introns_sort.bed` (131109,131110,131111)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.382GB.  
  PID: 131109;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 131111;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 131110;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_coverage.bed` (131120)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.382GB.  
  PID: 131120;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_coverage.bed` (131142)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.382GB.  
  PID: 131142;	Command: bedtools;	Return code: 0;	Memory used: 0.055GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/20.754829)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_rpkm.bed` (131168,131169,131170)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.382GB.  
  PID: 131168;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 131170;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 131169;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/20.754829)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_rpkm.bed` (131172,131173,131174)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 131172;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 131174;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 131173;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed` (131177,131178,131179)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 131177;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 131179;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 131178;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.23	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed --annotate` (131185)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 6.382GB.  
  PID: 131185;	Command: Rscript;	Return code: 0;	Memory used: 0.368GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_1_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_1_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed` (131206)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.382GB.  
  PID: 131206;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-27 12:31:54) elapsed: 70.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam` (131215)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.382GB.  
  PID: 131215;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (131222)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_1_plus_cuttrace_qd84090p'
Processing with 2 cores...
stdin is empty of data
Discarding 115 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr14_KI270724v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 80 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270438v1', 'chrUn_KI270515v1', 'chrUn_KI270539v1', 'chrUn_KI270587v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 80 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_exact_body_0-mer.bw'
Merging 80 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:23. Running peak memory: 6.382GB.  
  PID: 131222;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.542GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam` (132918)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.382GB.  
  PID: 132918;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (132927)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_1_minus_cuttrace_2ts1y31o'
Processing with 2 cores...
stdin is empty of data
stdin is empty of data
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_KI270722v1_random', 'chr14_KI270724v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:34. Running peak memory: 6.382GB.  
  PID: 132927;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.602GB

Starting cleanup: 77 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:25:37
*  Total elapsed time (all runs):  6:46:15
*         Peak memory (this run):  6.3822 GB
*        Pipeline completed time: 2020-02-27 12:51:09
