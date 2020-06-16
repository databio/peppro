### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_1 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-34c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/
*  Pipeline started at:   (06-15 07:17:16) elapsed: 1.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_1`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz

> `File_mb`	2469.73	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:17) elapsed: 2.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz` (261655)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 261655;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz` (261657)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 261657;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq` (261662)
<pre>
</pre>
Command completed. Elapsed time: 0:02:40. Running peak memory: 0.002GB.  
  PID: 261662;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2.fastq` (41970)
<pre>
</pre>
Command completed. Elapsed time: 0:01:34. Running peak memory: 0.002GB.  
  PID: 41970;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	97729860	PEPPRO	_RES_

> `Fastq_reads`	97729860	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz']

### FASTQ processing:  (06-15 07:23:33) elapsed: 376.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt` (233666)
<pre>
</pre>
Command completed. Elapsed time: 0:01:19. Running peak memory: 3.252GB.  
  PID: 233666;	Command: cutadapt;	Return code: 0;	Memory used: 3.252GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (233763,233764)
<pre>
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 3.252GB.  
  PID: 233763;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 233764;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	30221343	PEPPRO	_RES_

> `Trim_loss_rate`	69.08	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (234065)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-15 12:57:44) elapsed: 20051.0 _TIME_
Child process 234065 (fastqc) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/recover.lock.fastq__H9_PRO-seq_1_R1_processed.fastq
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/recover.lock.trimmed_fastqc

### Pipeline failed at:  (06-15 12:57:44) elapsed: 0.0 _TIME_

Total time: 5:40:29
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_1 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-12c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/
*  Pipeline started at:   (06-15 12:59:45) elapsed: 36.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_1`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz

> `File_mb`	2469.73	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 12:59:46) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz']

### FASTQ processing:  (06-15 12:59:46) elapsed: 0.0 _TIME_


> `cutadapt --version`
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/lock.fastq__H9_PRO-seq_1_R1_processed.fastq
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt` (14926)
<pre>
</pre>
Command completed. Elapsed time: 0:00:59. Running peak memory: 3.258GB.  
  PID: 14926;	Command: cutadapt;	Return code: 0;	Memory used: 3.258GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (15390,15391)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 3.258GB.  
  PID: 15390;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 15391;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	30221343	PEPPRO	_RES_

> `Trim_loss_rate`	69.08	PEPPRO	_RES_
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/lock.trimmed_fastqc
Overwriting target...
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (15447)
<pre>
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
Command completed. Elapsed time: 0:00:57. Running peak memory: 3.258GB.  
  PID: 15447;	Command: fastqc;	Return code: 0;	Memory used: 0.182GB

> `FastQC report r1`	fastqc/H9_PRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq` (15527)
<pre>
[INFO][0m 3516802 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 3.258GB.  
  PID: 15527;	Command: seqkit;	Return code: 0;	Memory used: 2.002GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq` (15609,15610)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.258GB.  
  PID: 15609;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 15610;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	34822522.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	17998987.0	PEPPRO	_RES_

> `Duplicate_reads`	3516802.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	36.8342	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/processed_R1.flag` (15713)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.258GB.  
  PID: 15713;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R2_cutadapt.txt` (15715)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 3.432GB.  
  PID: 15715;	Command: cutadapt;	Return code: 0;	Memory used: 3.432GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (15801,15802)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.432GB.  
  PID: 15801;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 15802;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	60442686	PEPPRO	_RES_

> `Trim_loss_rate`	38.15	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (16066)
<pre>
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
#
# A fatal error has been detected by the Java Runtime Environment:
#
#  SIGSEGV (0xb) at pc=0x00007f999e503631, pid=16066, tid=0x00007f9983bfb700
#
# JRE version: Java(TM) SE Runtime Environment (8.0_171-b11) (build 1.8.0_171-b11)
# Java VM: Java HotSpot(TM) 64-Bit Server VM (25.171-b11 mixed mode linux-amd64 compressed oops)
# Problematic frame:
# V  [libjvm.so+0x4b3631]  CompileBroker::invoke_compiler_on_method(CompileTask*)+0x5a1
#
# Core dump written. Default location: /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/core or core.16066
#
# An error report file with more information is saved as:
# /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc/hs_err_pid16066.log
#
# If you would like to submit a bug report, please visit:
#   http://bugreport.java.com/bugreport/crash.jsp
#
</pre>
Command completed. Elapsed time: 0:00:58. Running peak memory: 3.432GB.  
  PID: 16066;	Command: fastqc;	Return code: -6;	Memory used: 0.18GB

Subprocess returned nonzero result. Check above output for details
ERROR: Subprocess returned nonzero result, but pipeline is continuing because nofail=True
> `FastQC report r1`	fastqc/H9_PRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (16137)
<pre>
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
Command completed. Elapsed time: 0:00:54. Running peak memory: 3.432GB.  
  PID: 16137;	Command: fastqc;	Return code: 0;	Memory used: 0.17GB

> `FastQC report r2`	fastqc/H9_PRO-seq_1_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.histogram`  

> `fastq_pair -t 87956874 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq` (16204)
<pre>
Left paired: 30692275		Right paired: 30692275
Left single: 173668		Right single: 1557501
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:42. Running peak memory: 6.12GB.  
  PID: 16204;	Command: fastq_pair;	Return code: 0;	Memory used: 6.12GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq.paired.fq -o H9_PRO-seq_1 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt` (16448)
<pre>
</pre>
Command completed. Elapsed time: 0:00:59. Running peak memory: 6.12GB.  
  PID: 16448;	Command: flash;	Return code: 0;	Memory used: 0.094GB


### Plot adapter insertion distribution (06-15 13:10:55) elapsed: 669.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt -u 8` (16944)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.12GB.  
  PID: 16944;	Command: Rscript;	Return code: 0;	Memory used: 0.311GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 13:11:00) elapsed: 5.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `Degradation_ratio`	0.8896	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq` (16973)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 6.12GB.  
  PID: 16973;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/processed_R2.flag` (17341)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.12GB.  
  PID: 17341;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/repaired.flag`  

> `fastq_pair -t 87956874 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (17348)
<pre>
Left paired: 30003416		Right paired: 30003416
Left single: 217927		Right single: 1602455
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:32. Running peak memory: 6.178GB.  
  PID: 17348;	Command: fastq_pair;	Return code: 0;	Memory used: 6.178GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (17984)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 17984;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (17986)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 17986;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/repaired.flag` (17988)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 17988;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/dups_repaired.flag`  

> `fastq_pair -t 87956874 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq` (17989)
<pre>
Left paired: 26674489		Right paired: 26674489
Left single: 126830		Right single: 4931382
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:59. Running peak memory: 6.178GB.  
  PID: 17989;	Command: fastq_pair;	Return code: 0;	Memory used: 6.154GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq` (19059)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 19059;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq` (19060)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 19060;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/dups_repaired.flag` (19063)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 19063;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 13:15:55) elapsed: 295.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 13:15:55) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2` (19064)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 19064;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq` (19065)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
30003416 reads; of these:
  30003416 (100.00%) were unpaired; of these:
    27299076 (90.99%) aligned 0 times
    2704340 (9.01%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.01% overall alignment rate

> `Aligned_reads_human_rDNA`	5408680.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	8.95	PEPPRO	_RES_

### Map to human_rDNA (06-15 13:19:19) elapsed: 204.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2` (19503)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 19503;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R2.fq` (19504)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2704340 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 13:22:29) elapsed: 189.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_1 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/tmp909od977 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam` (20319,20320,20321)
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
[bam_sort_core] merging from 15 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:46:17. Running peak memory: 6.178GB.  
  PID: 20320;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20319;	Command: bowtie2;	Return code: 0;	Memory used: 3.802GB  
  PID: 20321;	Command: samtools;	Return code: 0;	Memory used: 0.897GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (24875)
<pre>
</pre>
Command completed. Elapsed time: 0:01:24. Running peak memory: 6.178GB.  
  PID: 24875;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	48339851	PEPPRO	_RES_

> `QC_filtered_reads`	27585022	PEPPRO	_RES_

> `Aligned_reads`	20754829.0	PEPPRO	_RES_

> `Alignment_rate`	34.34	PEPPRO	_RES_

> `Total_efficiency`	21.24	PEPPRO	_RES_

> `Read_depth`	3.26	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_1 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/tmp909od977 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam` (26249,26255,26256)
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
[bam_sort_core] merging from 14 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:41:14. Running peak memory: 6.178GB.  
  PID: 26249;	Command: bowtie2;	Return code: 0;	Memory used: 3.767GB  
  PID: 26255;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 26256;	Command: samtools;	Return code: 0;	Memory used: 0.897GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam` (30214)
<pre>
</pre>
Command completed. Elapsed time: 0:01:17. Running peak memory: 6.178GB.  
  PID: 30214;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 15:05:53) elapsed: 6204.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq` (30427)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 6.178GB.  
  PID: 30427;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq` (30459)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 6.178GB.  
  PID: 30459;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam` (30497)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 6.178GB.  
  PID: 30497;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	674956	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (30537)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.178GB.  
  PID: 30537;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/chr_sizes.bed` (30565,30566,30567,30568)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 30565;	Command: samtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 30567;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 30566;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 30568;	Command: grep;	Return code: 0;	Memory used: 0.003GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_noMT.bam` (30570)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 6.178GB.  
  PID: 30570;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (30622)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 30622;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (30623)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.178GB.  
  PID: 30623;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


### Split BAM file (06-15 15:08:57) elapsed: 184.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam` (30649,30650)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:46. Running peak memory: 6.178GB.  
  PID: 30649;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 30650;	Command: samtools;	Return code: 0;	Memory used: 5.83GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE2.bam` (31033,31034)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:27. Running peak memory: 6.178GB.  
  PID: 31034;	Command: samtools;	Return code: 0;	Memory used: 4.57GB  
  PID: 31033;	Command: samtools;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam` (31393)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 6.178GB.  
  PID: 31393;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam` (31436,31437)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:43. Running peak memory: 6.178GB.  
  PID: 31436;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 31437;	Command: samtools;	Return code: 0;	Memory used: 5.399GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE2.bam` (31655,31656)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:22. Running peak memory: 6.178GB.  
  PID: 31656;	Command: samtools;	Return code: 0;	Memory used: 4.278GB  
  PID: 31655;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


### Calculate library complexity (06-15 15:20:47) elapsed: 710.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam` (32034)
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
Command completed. Elapsed time: 0:01:41. Running peak memory: 6.178GB.  
  PID: 32034;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam` (32142)
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
.................................................................................._......._...........
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:47. Running peak memory: 6.178GB.  
  PID: 32142;	Command: preseq;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_counts.txt` (32237)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.178GB.  
  PID: 32237;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_plot` (32272)
<pre>
Processing H9_PRO-seq_1
INFO: Found real counts for H9_PRO-seq_1 - Total (M): 21.966708 Unique (M): 20.33023

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.178GB.  
  PID: 32272;	Command: Rscript;	Return code: 0;	Memory used: 0.211GB

> `Library complexity`	QC_hg38/H9_PRO-seq_1_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_1_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8774	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 15:24:50) elapsed: 243.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam` (32292)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 6.178GB.  
  PID: 32292;	Command: samtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv` (32424)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/tmp_H9_PRO-seq_1_PE1_b_xgcpnk'
Processing with 12 cores...
Discarding 96 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270720v1_random', 'chr14_KI270724v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 99 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270539v1', 'chrUn_KI270587v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.178GB.  
  PID: 32424;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.899GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	10983354.0	PEPPRO	_RES_

> `PBC2`	10983354.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_unmap.bam` (32480)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.178GB.  
  PID: 32480;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam`

> `Unmapped_reads`	6258301	PEPPRO	_RES_

### Split BAM by strand (06-15 15:25:48) elapsed: 58.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam` (32525)
<pre>
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 6.178GB.  
  PID: 32525;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam` (32589)
<pre>
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 6.178GB.  
  PID: 32589;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 15:28:07) elapsed: 139.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (32652)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 32652;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt` (32653)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.178GB.  
  PID: 32653;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.793GB


> `TSS_coding_score`	26.7	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt` (32686)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.178GB.  
  PID: 32686;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.925GB


> `TSS_non-coding_score`	10.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt` (32721)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.178GB.  
  PID: 32721;	Command: Rscript;	Return code: 0;	Memory used: 0.317GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_1_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_1_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt` (32742,32743,32744,32745)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 32742;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 32744;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 32743;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 32745;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt` (32747)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 32747;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-15 15:28:26) elapsed: 19.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed` (32750,32751)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.178GB.  
  PID: 32750;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 32751;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed` (32755,32756)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 32755;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 32756;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSS_density.bed` (32758,32759,32760,32761)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.178GB.  
  PID: 32759;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 32761;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 32758;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 32760;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_gene_body_density.bed` (32785,32786,32787)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.178GB.  
  PID: 32786;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 32785;	Command: bedtools;	Return code: 0;	Memory used: 0.063GB  
  PID: 32787;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/tmp6bum7qxi` (32822,32823,32824)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 32822;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 32824;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 32823;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/tmp6bum7qxi | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0210697) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/tmp6bum7qxi > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed` (32831)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 32831;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	26.45	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed` (32836)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.178GB.  
  PID: 32836;	Command: Rscript;	Return code: 0;	Memory used: 0.317GB

> `Pause index`	QC_hg38/H9_PRO-seq_1_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_1_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed` (32857)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 32857;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 15:29:32) elapsed: 65.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam`
20754829.0 7959339

> `Plus_FRiP`	0.38	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam`
20754829.0 7448642

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_gene_sort.bed` (33122,33123)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 33122;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 33123;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_gene_coverage.bed` (33126)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.178GB.  
  PID: 33126;	Command: bedtools;	Return code: 0;	Memory used: 0.069GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed.gz` (33158)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 33158;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed` (33159)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 33159;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 15:30:43) elapsed: 71.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed` (33168)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 33168;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed` (33170,33171,33172,33173)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 33170;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 33171;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 33173;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 33172;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed` (33175)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.178GB.  
  PID: 33175;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_minus_coverage.bed` (33193)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.178GB.  
  PID: 33193;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed` (33209,33210,33211,33212)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 33209;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 33210;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 33212;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 33211;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed` (33214)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.178GB.  
  PID: 33214;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_minus_coverage.bed` (33234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.178GB.  
  PID: 33234;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region"` (33251)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 33251;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed` (33252,33253,33254,33255)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 33252;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 33254;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 33253;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 33255;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (33257)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.178GB.  
  PID: 33257;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed` (33271)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.178GB.  
  PID: 33271;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR"` (33286)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 33286;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed` (33287,33288,33289,33290)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 33287;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 33288;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 33290;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 33289;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed` (33292)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.178GB.  
  PID: 33292;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_minus_coverage.bed` (33309)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.178GB.  
  PID: 33309;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR"` (33326)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 33326;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed` (33327,33328,33329,33330)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 33327;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 33328;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 33330;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 33329;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed` (33333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.178GB.  
  PID: 33333;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_minus_coverage.bed` (33349)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.178GB.  
  PID: 33349;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed` (33364,33365,33366,33367)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.178GB.  
  PID: 33364;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 33365;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 33367;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 33366;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed` (33373)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.178GB.  
  PID: 33373;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_minus_coverage.bed` (33390)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.178GB.  
  PID: 33390;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed` (33408,33409,33410,33411)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 33408;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 33410;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 33409;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 33411;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed` (33415)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.178GB.  
  PID: 33415;	Command: bedtools;	Return code: 0;	Memory used: 0.056GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_minus_coverage.bed` (33468)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.178GB.  
  PID: 33468;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB


### Plot cFRiF/FRiF (06-15 15:34:04) elapsed: 201.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_1 -z 3099922541 -n 11202246 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed` (33502)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 6.178GB.  
  PID: 33502;	Command: Rscript;	Return code: 0;	Memory used: 0.501GB

> `cFRiF`	QC_hg38/H9_PRO-seq_1_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_1_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_1 -z 3099922541 -n 11202246 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed` (33548)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.178GB.  
  PID: 33548;	Command: Rscript;	Return code: 0;	Memory used: 0.501GB

> `FRiF`	QC_hg38/H9_PRO-seq_1_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_1_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 15:35:02) elapsed: 58.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_exons_sort.bed` (33698,33699)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.178GB.  
  PID: 33698;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 33699;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_introns_sort.bed` (33705,33706,33707)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.178GB.  
  PID: 33705;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 33707;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 33706;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_coverage.bed` (33719)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.178GB.  
  PID: 33719;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_coverage.bed` (33747)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.178GB.  
  PID: 33747;	Command: bedtools;	Return code: 0;	Memory used: 0.055GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/20.754829)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_rpkm.bed` (33772,33773,33774)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 33772;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 33774;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 33773;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/20.754829)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_rpkm.bed` (33776,33777,33778)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.178GB.  
  PID: 33776;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 33778;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 33777;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed` (33781,33782,33783)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 33781;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 33783;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 33782;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.23	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed --annotate` (33789)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.178GB.  
  PID: 33789;	Command: Rscript;	Return code: 0;	Memory used: 0.317GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_1_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_1_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed` (33810)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.178GB.  
  PID: 33810;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 15:36:13) elapsed: 71.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam` (33818)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 6.178GB.  
  PID: 33818;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 20754829.0` (33829)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_1_plus_cuttrace_exvi9dtb'
Processing with 4 cores...
stdin is empty of data
Discarding 115 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr14_KI270724v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 80 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270438v1', 'chrUn_KI270515v1', 'chrUn_KI270539v1', 'chrUn_KI270587v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 80 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_exact_body_0-mer.bw'
Merging 80 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:18. Running peak memory: 6.178GB.  
  PID: 33829;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.346GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam` (35275)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.178GB.  
  PID: 35275;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 20754829.0` (35283)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_1_minus_cuttrace_j2gcnfml'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_KI270722v1_random', 'chr14_KI270724v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:12. Running peak memory: 6.178GB.  
  PID: 35283;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.763GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:50:51
*  Total elapsed time (all runs):  5:20:58
*         Peak memory (this run):  6.178 GB
*        Pipeline completed time: 2020-06-15 15:50:01
