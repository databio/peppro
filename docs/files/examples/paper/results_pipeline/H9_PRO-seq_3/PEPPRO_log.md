### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_3 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-34c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/
*  Pipeline started at:   (06-15 07:17:16) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_3`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz

> `File_mb`	2340.86	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:17) elapsed: 2.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz` (261654)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 261654;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz` (261658)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 261658;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq` (261663)
<pre>
</pre>
Command completed. Elapsed time: 0:02:18. Running peak memory: 0.002GB.  
  PID: 261663;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2.fastq` (27124)
<pre>
</pre>
Command completed. Elapsed time: 0:02:15. Running peak memory: 0.002GB.  
  PID: 27124;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	92628962	PEPPRO	_RES_

> `Fastq_reads`	92628962	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz']

### FASTQ processing:  (06-15 07:23:30) elapsed: 373.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt` (233633)
<pre>
</pre>
Command completed. Elapsed time: 0:01:31. Running peak memory: 3.239GB.  
  PID: 233633;	Command: cutadapt;	Return code: 0;	Memory used: 3.239GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (233962,233963)
<pre>
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 3.239GB.  
  PID: 233962;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 233963;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	27419170	PEPPRO	_RES_

> `Trim_loss_rate`	70.4	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (234067)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-15 12:57:37) elapsed: 20046.0 _TIME_
Child process 234067 (fastqc) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/recover.lock.fastq__H9_PRO-seq_3_R1_processed.fastq
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/recover.lock.trimmed_fastqc

### Pipeline failed at:  (06-15 12:57:37) elapsed: 0.0 _TIME_

Total time: 5:40:22
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_3 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-33c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/
*  Pipeline started at:   (06-15 12:59:45) elapsed: 35.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_3`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz

> `File_mb`	2340.86	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 12:59:46) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz']

### FASTQ processing:  (06-15 12:59:46) elapsed: 0.0 _TIME_


> `cutadapt --version`
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/lock.fastq__H9_PRO-seq_3_R1_processed.fastq
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt` (374018)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 3.233GB.  
  PID: 374018;	Command: cutadapt;	Return code: 0;	Memory used: 3.233GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (374356,374357)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.233GB.  
  PID: 374356;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 374357;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	27419170	PEPPRO	_RES_

> `Trim_loss_rate`	70.4	PEPPRO	_RES_
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/lock.trimmed_fastqc
Overwriting target...
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (374400)
<pre>
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
Command completed. Elapsed time: 0:01:00. Running peak memory: 3.233GB.  
  PID: 374400;	Command: fastqc;	Return code: 0;	Memory used: 0.181GB

> `FastQC report r1`	fastqc/H9_PRO-seq_3_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq` (374469)
<pre>
[INFO][0m 1862033 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 3.233GB.  
  PID: 374469;	Command: seqkit;	Return code: 0;	Memory used: 2.006GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq` (374528,374529)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.233GB.  
  PID: 374528;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 374529;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	34725419.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	18106301.0	PEPPRO	_RES_

> `Duplicate_reads`	1862033.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	39.0943	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/processed_R1.flag` (374598)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.233GB.  
  PID: 374598;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R2_cutadapt.txt` (374600)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 3.463GB.  
  PID: 374600;	Command: cutadapt;	Return code: 0;	Memory used: 3.463GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (374674,374676)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.463GB.  
  PID: 374674;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 374676;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	54838340	PEPPRO	_RES_

> `Trim_loss_rate`	40.8	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (374900)
<pre>
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
Command completed. Elapsed time: 0:01:02. Running peak memory: 3.463GB.  
  PID: 374900;	Command: fastqc;	Return code: 0;	Memory used: 0.178GB

> `FastQC report r1`	fastqc/H9_PRO-seq_3_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (374971)
<pre>
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
Command completed. Elapsed time: 0:00:57. Running peak memory: 3.463GB.  
  PID: 374971;	Command: fastqc;	Return code: 0;	Memory used: 0.171GB

> `FastQC report r2`	fastqc/H9_PRO-seq_3_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.histogram`  

> `fastq_pair -t 83366065 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq` (375045)
<pre>
Left paired: 28081513		Right paired: 28081513
Left single: 126667		Right single: 1382160
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:58. Running peak memory: 6.109GB.  
  PID: 375045;	Command: fastq_pair;	Return code: 0;	Memory used: 6.109GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq.paired.fq -o H9_PRO-seq_3 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt` (375234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 6.109GB.  
  PID: 375234;	Command: flash;	Return code: 0;	Memory used: 0.098GB


### Plot adapter insertion distribution (06-15 13:10:57) elapsed: 671.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt -u 8` (375678)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.109GB.  
  PID: 375678;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_3_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_3_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 13:11:02) elapsed: 5.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `Degradation_ratio`	0.9769	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq` (375707)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 6.109GB.  
  PID: 375707;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/processed_R2.flag` (375726)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 375726;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/repaired.flag`  

> `fastq_pair -t 83366065 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (375727)
<pre>
Left paired: 27272197		Right paired: 27272197
Left single: 146973		Right single: 1431712
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:03:15. Running peak memory: 6.109GB.  
  PID: 375727;	Command: fastq_pair;	Return code: 0;	Memory used: 5.962GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (375897)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 375897;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (375898)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 375898;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/repaired.flag` (375900)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 375900;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/dups_repaired.flag`  

> `fastq_pair -t 83366065 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq` (375901)
<pre>
Left paired: 25551175		Right paired: 25551175
Left single: 94450		Right single: 3152734
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:45. Running peak memory: 6.109GB.  
  PID: 375901;	Command: fastq_pair;	Return code: 0;	Memory used: 5.253GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq` (376242)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 376242;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq` (376244)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 376244;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/dups_repaired.flag` (376245)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 376245;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 13:17:26) elapsed: 384.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 13:17:26) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2` (376246)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 376246;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq` (376247)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_3 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
27272197 reads; of these:
  27272197 (100.00%) were unpaired; of these:
    24569260 (90.09%) aligned 0 times
    2702937 (9.91%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.91% overall alignment rate

> `Aligned_reads_human_rDNA`	5405874.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.86	PEPPRO	_RES_

### Map to human_rDNA (06-15 13:20:55) elapsed: 208.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2` (377186)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 377186;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R2.fq` (377187)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_3 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2702937 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 13:24:53) elapsed: 238.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_3 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/tmpofh5asc7 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam` (381212,381213,381214)
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
[bam_sort_core] merging from 14 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:50:45. Running peak memory: 6.109GB.  
  PID: 381213;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 381212;	Command: bowtie2;	Return code: 0;	Memory used: 3.765GB  
  PID: 381214;	Command: samtools;	Return code: 0;	Memory used: 0.898GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (395318)
<pre>
</pre>
Command completed. Elapsed time: 0:02:09. Running peak memory: 6.109GB.  
  PID: 395318;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	43295089	PEPPRO	_RES_

> `QC_filtered_reads`	25127415	PEPPRO	_RES_

> `Aligned_reads`	18167674.5	PEPPRO	_RES_

> `Alignment_rate`	33.13	PEPPRO	_RES_

> `Total_efficiency`	19.61	PEPPRO	_RES_

> `Read_depth`	2.98	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_3 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/tmpofh5asc7 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam` (401719,401724,401725)
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
[bam_sort_core] merging from 13 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:47:32. Running peak memory: 6.109GB.  
  PID: 401719;	Command: bowtie2;	Return code: 0;	Memory used: 3.762GB  
  PID: 401724;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 401725;	Command: samtools;	Return code: 0;	Memory used: 0.921GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam` (454442)
<pre>
</pre>
Command completed. Elapsed time: 0:01:26. Running peak memory: 6.109GB.  
  PID: 454442;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 15:20:43) elapsed: 6950.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq` (21123)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 6.109GB.  
  PID: 21123;	Command: pigz;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq` (28309)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.109GB.  
  PID: 28309;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam` (30700)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 6.109GB.  
  PID: 30700;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	711009	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (35171)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.109GB.  
  PID: 35171;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/chr_sizes.bed` (35239,35240,35241,35242)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 35241;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 35239;	Command: samtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 35242;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 35240;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_noMT.bam` (35244)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 6.109GB.  
  PID: 35244;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (35321)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 35321;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (35322)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.109GB.  
  PID: 35322;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


### Split BAM file (06-15 15:23:45) elapsed: 182.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam` (38421,38436)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:41. Running peak memory: 6.109GB.  
  PID: 38421;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 38436;	Command: samtools;	Return code: 0;	Memory used: 5.156GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE2.bam` (77855,77860)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:17. Running peak memory: 6.109GB.  
  PID: 77860;	Command: samtools;	Return code: 0;	Memory used: 4.019GB  
  PID: 77855;	Command: samtools;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam` (79965)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 6.109GB.  
  PID: 79965;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam` (80272,80273)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:42. Running peak memory: 6.109GB.  
  PID: 80272;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 80273;	Command: samtools;	Return code: 0;	Memory used: 5.033GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE2.bam` (118983,118984)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:10. Running peak memory: 6.109GB.  
  PID: 118983;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 118984;	Command: samtools;	Return code: 0;	Memory used: 4.074GB


### Calculate library complexity (06-15 15:35:16) elapsed: 690.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam` (119394)
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
Command completed. Elapsed time: 0:01:47. Running peak memory: 6.109GB.  
  PID: 119394;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam` (157921)
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
.............................................._.........._............................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:02:04. Running peak memory: 6.109GB.  
  PID: 157921;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_counts.txt` (191204)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.109GB.  
  PID: 191204;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_plot` (191240)
<pre>
Processing H9_PRO-seq_3
INFO: Found real counts for H9_PRO-seq_3 - Total (M): 19.481041 Unique (M): 19.005606

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.109GB.  
  PID: 191240;	Command: Rscript;	Return code: 0;	Memory used: 0.273GB

> `Library complexity`	QC_hg38/H9_PRO-seq_3_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_3_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8821	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 15:39:44) elapsed: 269.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam` (191263)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 6.109GB.  
  PID: 191263;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv` (191278)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/tmp_H9_PRO-seq_3_PE1_svpmvopz'
Processing with 12 cores...
Discarding 101 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 94 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 6.109GB.  
  PID: 191278;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.51GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	9740520.5	PEPPRO	_RES_

> `PBC2`	9740520.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_unmap.bam` (191640)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.109GB.  
  PID: 191640;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam`

> `Unmapped_reads`	5843431	PEPPRO	_RES_

### Split BAM by strand (06-15 15:40:41) elapsed: 57.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam` (191711)
<pre>
</pre>
Command completed. Elapsed time: 0:01:19. Running peak memory: 6.109GB.  
  PID: 191711;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam` (191798)
<pre>
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 6.109GB.  
  PID: 191798;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 15:43:12) elapsed: 151.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (201242)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 201242;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt` (201333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.109GB.  
  PID: 201333;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.723GB


> `TSS_coding_score`	28.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt` (205245)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.109GB.  
  PID: 205245;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.763GB


> `TSS_non-coding_score`	10.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt` (209694)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.109GB.  
  PID: 209694;	Command: Rscript;	Return code: 0;	Memory used: 0.123GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_3_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_3_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt` (220364,220380,220403,220416)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 220364;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 220403;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 220380;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 220416;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt` (220478)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 220478;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-15 15:43:39) elapsed: 27.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed` (220576,220582)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.109GB.  
  PID: 220576;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 220582;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed` (221454,221455)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 221454;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 221455;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSS_density.bed` (221577,221587,221590,221596)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.109GB.  
  PID: 221587;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 221596;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 221577;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 221590;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_gene_body_density.bed` (228579,228581,228582)
<pre>
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 6.109GB.  
  PID: 228581;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 228579;	Command: bedtools;	Return code: 0;	Memory used: 0.071GB  
  PID: 228582;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/tmpg2fp3lds` (230727,230728,230729)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 230727;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 230729;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 230728;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/tmpg2fp3lds | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0188796) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/tmpg2fp3lds > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed` (230735)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 230735;	Command: awk;	Return code: 0;	Memory used: 0.004GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	27.32	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed` (230740)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.109GB.  
  PID: 230740;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Pause index`	QC_hg38/H9_PRO-seq_3_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_3_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed` (230781)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 230781;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 15:45:31) elapsed: 112.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam`
18167674.5 7066591

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam`
18167674.5 6601482

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_gene_sort.bed` (230862,230863)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 230862;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 230863;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_gene_coverage.bed` (230866)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 6.109GB.  
  PID: 230866;	Command: bedtools;	Return code: 0;	Memory used: 0.063GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed.gz` (230902)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 230902;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed` (230903)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 230903;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 15:46:47) elapsed: 76.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed` (230911)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 230911;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed` (230913,230914,230915,230916)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 230913;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 230914;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 230916;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 230915;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed` (230919)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.109GB.  
  PID: 230919;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_minus_coverage.bed` (230932)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.109GB.  
  PID: 230932;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed` (230944,230945,230946,230947)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 230944;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 230945;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 230947;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 230946;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed` (230950)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.109GB.  
  PID: 230950;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_minus_coverage.bed` (245408)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.109GB.  
  PID: 245408;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region"` (264353)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 264353;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed` (264409,264414,264421,264426)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 264409;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 264421;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 264414;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 264426;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed` (265975)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.109GB.  
  PID: 265975;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_minus_coverage.bed` (269452)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.109GB.  
  PID: 269452;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR"` (269467)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 269467;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed` (269468,269470,269471,269472)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 269468;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 269470;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 269472;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 269471;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed` (269474)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.109GB.  
  PID: 269474;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_minus_coverage.bed` (269488)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.109GB.  
  PID: 269488;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR"` (269516)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 269516;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed` (269517,269518,269519,269520)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 269517;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 269518;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 269520;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB  
  PID: 269519;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed` (269523)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.109GB.  
  PID: 269523;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_minus_coverage.bed` (269537)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.109GB.  
  PID: 269537;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed` (269562,269563,269564,269565)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.109GB.  
  PID: 269562;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 269563;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 269565;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 269564;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed` (269570)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 6.109GB.  
  PID: 269570;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_minus_coverage.bed` (269751)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.109GB.  
  PID: 269751;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed` (282124,282249,282254,282255)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.109GB.  
  PID: 282124;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 282254;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 282249;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 282255;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed` (283800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.109GB.  
  PID: 283800;	Command: bedtools;	Return code: 0;	Memory used: 0.058GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_minus_coverage.bed` (299691)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 6.109GB.  
  PID: 299691;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB


### Plot cFRiF/FRiF (06-15 15:50:43) elapsed: 236.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_3 -z 3099922541 -n 9939223 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed` (308434)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 6.109GB.  
  PID: 308434;	Command: Rscript;	Return code: 0;	Memory used: 0.557GB

> `cFRiF`	QC_hg38/H9_PRO-seq_3_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_3_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_3 -z 3099922541 -n 9939223 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed` (308508)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 6.109GB.  
  PID: 308508;	Command: Rscript;	Return code: 0;	Memory used: 0.442GB

> `FRiF`	QC_hg38/H9_PRO-seq_3_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_3_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 15:52:27) elapsed: 104.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_exons_sort.bed` (308606,308607)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.109GB.  
  PID: 308607;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB  
  PID: 308606;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_introns_sort.bed` (308613,308614,308615)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.109GB.  
  PID: 308613;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 308615;	Command: bedtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 308614;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_coverage.bed` (308622)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.109GB.  
  PID: 308622;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_coverage.bed` (308647)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 6.109GB.  
  PID: 308647;	Command: bedtools;	Return code: 0;	Memory used: 0.031GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.1676745)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_rpkm.bed` (316390,316391,316392)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 316390;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 316392;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 316391;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.1676745)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_rpkm.bed` (316783,316796,316797)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.109GB.  
  PID: 316783;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 316797;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 316796;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed` (317244,317260,317268)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 317244;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 317268;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 317260;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.28	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed --annotate` (317508)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 6.109GB.  
  PID: 317508;	Command: Rscript;	Return code: 0;	Memory used: 0.11GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_3_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_3_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed` (322436)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.109GB.  
  PID: 322436;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 15:54:00) elapsed: 93.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam` (322534)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.109GB.  
  PID: 322534;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 18167674.5` (327192)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_3_plus_cuttrace_wm7hae5s'
Processing with 4 cores...
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:08:15. Running peak memory: 6.109GB.  
  PID: 327192;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.595GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam` (415919)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.109GB.  
  PID: 415919;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 18167674.5` (415926)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_3_minus_cuttrace_bf2ykl7_'
Processing with 4 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_exact_body_0-mer.bw'
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:48. Running peak memory: 6.109GB.  
  PID: 415926;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.752GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:11:09
*  Total elapsed time (all runs):  5:59:47
*         Peak memory (this run):  6.1094 GB
*        Pipeline completed time: 2020-06-15 16:10:19
