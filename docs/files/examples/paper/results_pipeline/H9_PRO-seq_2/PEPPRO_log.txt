### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_2 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-34c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/
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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_2`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz

> `File_mb`	2754.18	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:17) elapsed: 2.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz` (261653)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 261653;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz` (261656)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 261656;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq` (261659)
<pre>
</pre>
Command completed. Elapsed time: 0:03:01. Running peak memory: 0.002GB.  
  PID: 261659;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2.fastq` (176590)
<pre>
</pre>
Command completed. Elapsed time: 0:02:06. Running peak memory: 0.002GB.  
  PID: 176590;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	115714026	PEPPRO	_RES_

> `Fastq_reads`	115714026	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz']

### FASTQ processing:  (06-15 07:25:17) elapsed: 479.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt` (233983)
<pre>
</pre>
Command completed. Elapsed time: 0:01:59. Running peak memory: 3.774GB.  
  PID: 233983;	Command: cutadapt;	Return code: 0;	Memory used: 3.774GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (234120,234121)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 3.774GB.  
  PID: 234120;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 234121;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	28329255	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (234153)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-15 12:58:05) elapsed: 19969.0 _TIME_
Child process 234153 (fastqc) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/recover.lock.fastq__H9_PRO-seq_2_R1_processed.fastq
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/recover.lock.trimmed_fastqc

### Pipeline failed at:  (06-15 12:58:05) elapsed: 0.0 _TIME_

Total time: 5:40:51
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_2 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-28c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/
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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_2`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz

> `File_mb`	2754.18	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 12:59:46) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz']

### FASTQ processing:  (06-15 12:59:46) elapsed: 0.0 _TIME_


> `cutadapt --version`
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/lock.fastq__H9_PRO-seq_2_R1_processed.fastq
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt` (190216)
<pre>
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 3.915GB.  
  PID: 190216;	Command: cutadapt;	Return code: 0;	Memory used: 3.915GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (191521,191522)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.915GB.  
  PID: 191521;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 191522;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	28329255	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/lock.trimmed_fastqc
Overwriting target...
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (191586)
<pre>
Started analysis of H9_PRO-seq_2_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_2_R1_processed.fastq
Analysis complete for H9_PRO-seq_2_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 3.915GB.  
  PID: 191586;	Command: fastqc;	Return code: 0;	Memory used: 0.18GB

> `FastQC report r1`	fastqc/H9_PRO-seq_2_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq` (192172)
<pre>
[INFO][0m 3047417 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 3.915GB.  
  PID: 192172;	Command: seqkit;	Return code: 0;	Memory used: 1.03GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq` (192734,192735)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.915GB.  
  PID: 192734;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 192735;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	45847892.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	28610139.0	PEPPRO	_RES_

> `Duplicate_reads`	3047417.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.4497	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/processed_R1.flag` (193051)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.915GB.  
  PID: 193051;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R2_cutadapt.txt` (193053)
<pre>
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 3.915GB.  
  PID: 193053;	Command: cutadapt;	Return code: 0;	Memory used: 3.373GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (193811,193812)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.915GB.  
  PID: 193811;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 193812;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	56658510	PEPPRO	_RES_

> `Trim_loss_rate`	51.04	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (194145)
<pre>
Started analysis of H9_PRO-seq_2_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_2_R1_processed.fastq
Analysis complete for H9_PRO-seq_2_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:56. Running peak memory: 3.915GB.  
  PID: 194145;	Command: fastqc;	Return code: 0;	Memory used: 0.178GB

> `FastQC report r1`	fastqc/H9_PRO-seq_2_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (194953)
<pre>
Started analysis of H9_PRO-seq_2_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_2_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_2_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 3.915GB.  
  PID: 194953;	Command: fastqc;	Return code: 0;	Memory used: 0.17GB

> `FastQC report r2`	fastqc/H9_PRO-seq_2_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.histogram`  

> `fastq_pair -t 104142623 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq` (195894)
<pre>
Left paired: 29044502		Right paired: 29044502
Left single: 202372		Right single: 2678154
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:21. Running peak memory: 6.128GB.  
  PID: 195894;	Command: fastq_pair;	Return code: 0;	Memory used: 6.128GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq.paired.fq -o H9_PRO-seq_2 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt` (198360)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 6.128GB.  
  PID: 198360;	Command: flash;	Return code: 0;	Memory used: 0.09GB


### Plot adapter insertion distribution (06-15 13:10:01) elapsed: 615.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt -u 8` (199407)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.128GB.  
  PID: 199407;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_2_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_2_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 13:10:10) elapsed: 9.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `Degradation_ratio`	1.0321	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq` (199588)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.128GB.  
  PID: 199588;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/processed_R2.flag` (201244)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.128GB.  
  PID: 201244;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/repaired.flag`  

> `fastq_pair -t 104142623 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (201245)
<pre>
Left paired: 28105674		Right paired: 28105674
Left single: 223581		Right single: 2750162
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:14. Running peak memory: 6.281GB.  
  PID: 201245;	Command: fastq_pair;	Return code: 0;	Memory used: 6.281GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (203149)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 203149;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (203154)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 203154;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/repaired.flag` (203155)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 203155;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/dups_repaired.flag`  

> `fastq_pair -t 104142623 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq` (203156)
<pre>
Left paired: 25286654		Right paired: 25286654
Left single: 163446		Right single: 5569182
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:48. Running peak memory: 6.281GB.  
  PID: 203156;	Command: fastq_pair;	Return code: 0;	Memory used: 5.417GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq` (209107)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 209107;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq` (209132)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 209132;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/dups_repaired.flag` (209133)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 209133;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 13:14:30) elapsed: 261.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 13:14:30) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2` (209137)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 209137;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq` (209138)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_2 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
28105674 reads; of these:
  28105674 (100.00%) were unpaired; of these:
    25192086 (89.63%) aligned 0 times
    2913588 (10.37%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.37% overall alignment rate

> `Aligned_reads_human_rDNA`	5827176.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.28	PEPPRO	_RES_

### Map to human_rDNA (06-15 13:17:14) elapsed: 164.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2` (214548)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 214548;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R2.fq` (214549)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_2 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2913588 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 13:19:45) elapsed: 150.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_2 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/tmpy4mio3xb -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam` (216935,216936,216937)
<pre>
2252266 reads skipped
0 reads lost
25192086 reads; of these:
  25192086 (100.00%) were paired; of these:
    10443566 (41.46%) aligned concordantly 0 times
    12095912 (48.01%) aligned concordantly exactly 1 time
    2652608 (10.53%) aligned concordantly >1 times
    ----
    10443566 pairs aligned concordantly 0 times; of these:
      2323449 (22.25%) aligned discordantly 1 time
    ----
    8120117 pairs aligned 0 times concordantly or discordantly; of these:
      16240234 mates make up the pairs; of these:
        6664276 (41.04%) aligned 0 times
        3498753 (21.54%) aligned exactly 1 time
        6077205 (37.42%) aligned >1 times
86.77% overall alignment rate
[bam_sort_core] merging from 14 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:41:27. Running peak memory: 6.281GB.  
  PID: 216935;	Command: bowtie2;	Return code: 0;	Memory used: 3.778GB  
  PID: 216936;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 216937;	Command: samtools;	Return code: 0;	Memory used: 0.924GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (222569)
<pre>
</pre>
Command completed. Elapsed time: 0:01:25. Running peak memory: 6.281GB.  
  PID: 222569;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	43719896	PEPPRO	_RES_

> `QC_filtered_reads`	25554169	PEPPRO	_RES_

> `Aligned_reads`	18165727.0	PEPPRO	_RES_

> `Alignment_rate`	32.06	PEPPRO	_RES_

> `Total_efficiency`	15.7	PEPPRO	_RES_

> `Read_depth`	3.11	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_2 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/tmpy4mio3xb -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam` (223921,223934,223935)
<pre>
23034388 reads; of these:
  23034388 (100.00%) were paired; of these:
    9342849 (40.56%) aligned concordantly 0 times
    11259901 (48.88%) aligned concordantly exactly 1 time
    2431638 (10.56%) aligned concordantly >1 times
    ----
    9342849 pairs aligned concordantly 0 times; of these:
      2170368 (23.23%) aligned discordantly 1 time
    ----
    7172481 pairs aligned 0 times concordantly or discordantly; of these:
      14344962 mates make up the pairs; of these:
        6003927 (41.85%) aligned 0 times
        3261117 (22.73%) aligned exactly 1 time
        5079918 (35.41%) aligned >1 times
86.97% overall alignment rate
[bam_sort_core] merging from 13 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:37:34. Running peak memory: 6.281GB.  
  PID: 223934;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 223921;	Command: bowtie2;	Return code: 0;	Memory used: 3.763GB  
  PID: 223935;	Command: samtools;	Return code: 0;	Memory used: 0.906GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam` (227336)
<pre>
</pre>
Command completed. Elapsed time: 0:01:16. Running peak memory: 6.281GB.  
  PID: 227336;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 14:53:11) elapsed: 5606.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq` (227422)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 6.281GB.  
  PID: 227422;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq` (227458)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 6.281GB.  
  PID: 227458;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam` (227500)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 6.281GB.  
  PID: 227500;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	840435	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (227533)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.281GB.  
  PID: 227533;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/chr_sizes.bed` (227564,227565,227566,227567)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 227566;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 227564;	Command: samtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 227567;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 227565;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_noMT.bam` (227569)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.281GB.  
  PID: 227569;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (227733)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 227733;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (227734)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 6.281GB.  
  PID: 227734;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


### Split BAM file (06-15 14:55:54) elapsed: 163.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam` (227759,227760)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:25. Running peak memory: 6.281GB.  
  PID: 227760;	Command: samtools;	Return code: 0;	Memory used: 4.811GB  
  PID: 227759;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE2.bam` (227924,227925)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:59. Running peak memory: 6.281GB.  
  PID: 227924;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 227925;	Command: samtools;	Return code: 0;	Memory used: 4.146GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam` (228332)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.281GB.  
  PID: 228332;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam` (228365,228366)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:23. Running peak memory: 6.281GB.  
  PID: 228366;	Command: samtools;	Return code: 0;	Memory used: 4.628GB  
  PID: 228365;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE2.bam` (228533,228534)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:57. Running peak memory: 6.281GB.  
  PID: 228533;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 228534;	Command: samtools;	Return code: 0;	Memory used: 3.933GB


### Calculate library complexity (06-15 15:05:59) elapsed: 605.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam` (228831)
<pre>
BAM_INPUT
TOTAL READS     = 18471612
COUNTS_SUM      = 18471612
DISTINCT READS  = 1.4838e+07
DISTINCT COUNTS = 324
MAX COUNT       = 27760
COUNTS OF 1     = 1.31407e+07
OBSERVED COUNTS (27761)
1	13140744
2	1146563
3	271114
4	106759
5	54977
6	32098
7	20377
8	13755
9	9695
10	7243
11	5405
12	4178
13	3366
14	2630
15	2184
16	1670
17	1556
18	1258
19	1062
20	922
21	848
22	719
23	630
24	559
25	515
26	451
27	428
28	367
29	288
30	299
31	270
32	267
33	217
34	243
35	187
36	184
37	185
38	182
39	115
40	118
41	128
42	138
43	109
44	108
45	116
46	94
47	93
48	89
49	86
50	74
51	68
52	83
53	92
54	70
55	66
56	62
57	55
58	57
59	51
60	52
61	61
62	42
63	45
64	39
65	35
66	29
67	50
68	36
69	25
70	39
71	35
72	28
73	31
74	32
75	27
76	27
77	24
78	27
79	28
80	21
81	23
82	16
83	21
84	29
85	23
86	16
87	19
88	19
89	17
90	12
91	10
92	17
93	11
94	17
95	10
96	17
97	22
98	19
99	17
100	9
101	14
102	19
103	10
104	8
105	9
106	9
107	16
108	8
109	8
110	5
111	8
112	5
113	13
114	12
115	8
116	12
117	7
118	4
119	9
120	7
121	6
122	5
123	11
124	8
125	2
126	5
127	5
128	7
129	6
130	4
131	5
132	6
133	5
134	4
135	2
136	3
137	4
138	11
139	5
140	9
141	3
142	2
143	5
144	5
145	4
146	2
147	6
148	2
149	5
150	5
151	2
152	4
153	3
154	3
155	4
156	5
157	4
158	5
159	1
160	3
161	4
162	6
163	5
164	2
165	5
166	6
167	1
168	1
169	2
170	2
171	3
172	3
173	1
174	3
175	2
176	1
177	5
178	2
179	3
180	2
181	4
183	1
184	2
186	3
187	1
188	1
190	3
192	4
194	4
196	3
197	2
198	3
199	2
200	3
201	2
202	1
203	1
204	4
206	1
207	4
208	4
210	3
212	3
213	3
215	2
219	1
220	3
221	2
222	2
223	1
224	2
225	3
226	3
228	1
230	1
232	1
233	2
234	1
235	2
236	2
238	2
240	1
242	2
244	2
246	1
249	1
253	2
255	1
256	1
257	1
258	2
259	2
261	1
263	2
265	1
268	1
269	2
270	1
273	1
277	1
279	1
282	1
283	1
284	1
290	1
292	1
299	1
302	2
303	1
306	1
307	1
310	1
311	1
312	1
314	1
315	1
317	1
318	1
319	2
320	1
321	1
322	1
325	1
328	2
331	1
335	2
338	1
341	1
345	1
350	1
355	1
360	2
370	1
373	1
374	1
390	2
391	1
394	1
397	1
412	1
415	1
426	1
429	1
442	1
444	1
452	1
488	1
494	2
502	1
505	1
509	1
515	1
549	2
551	1
584	1
694	1
701	1
703	1
708	1
710	1
712	1
714	1
809	1
881	1
918	1
962	1
1001	1
1158	1
1264	1
1324	1
1628	1
1848	1
1890	1
1899	1
2227	1
2457	1
2584	1
3405	1
4384	1
9005	1
10491	1
13480	1
23174	1
27760	1

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
</pre>
Command completed. Elapsed time: 0:01:33. Running peak memory: 6.281GB.  
  PID: 228831;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam` (228910)
<pre>
BAM_INPUT
TOTAL READS     = 18471612
DISTINCT READS  = 1.4838e+07
DISTINCT COUNTS = 324
MAX COUNT       = 27760
COUNTS OF 1     = 1.31407e+07
MAX TERMS       = 100
OBSERVED COUNTS (27761)
1	13140744
2	1146563
3	271114
4	106759
5	54977
6	32098
7	20377
8	13755
9	9695
10	7243
11	5405
12	4178
13	3366
14	2630
15	2184
16	1670
17	1556
18	1258
19	1062
20	922
21	848
22	719
23	630
24	559
25	515
26	451
27	428
28	367
29	288
30	299
31	270
32	267
33	217
34	243
35	187
36	184
37	185
38	182
39	115
40	118
41	128
42	138
43	109
44	108
45	116
46	94
47	93
48	89
49	86
50	74
51	68
52	83
53	92
54	70
55	66
56	62
57	55
58	57
59	51
60	52
61	61
62	42
63	45
64	39
65	35
66	29
67	50
68	36
69	25
70	39
71	35
72	28
73	31
74	32
75	27
76	27
77	24
78	27
79	28
80	21
81	23
82	16
83	21
84	29
85	23
86	16
87	19
88	19
89	17
90	12
91	10
92	17
93	11
94	17
95	10
96	17
97	22
98	19
99	17
100	9
101	14
102	19
103	10
104	8
105	9
106	9
107	16
108	8
109	8
110	5
111	8
112	5
113	13
114	12
115	8
116	12
117	7
118	4
119	9
120	7
121	6
122	5
123	11
124	8
125	2
126	5
127	5
128	7
129	6
130	4
131	5
132	6
133	5
134	4
135	2
136	3
137	4
138	11
139	5
140	9
141	3
142	2
143	5
144	5
145	4
146	2
147	6
148	2
149	5
150	5
151	2
152	4
153	3
154	3
155	4
156	5
157	4
158	5
159	1
160	3
161	4
162	6
163	5
164	2
165	5
166	6
167	1
168	1
169	2
170	2
171	3
172	3
173	1
174	3
175	2
176	1
177	5
178	2
179	3
180	2
181	4
183	1
184	2
186	3
187	1
188	1
190	3
192	4
194	4
196	3
197	2
198	3
199	2
200	3
201	2
202	1
203	1
204	4
206	1
207	4
208	4
210	3
212	3
213	3
215	2
219	1
220	3
221	2
222	2
223	1
224	2
225	3
226	3
228	1
230	1
232	1
233	2
234	1
235	2
236	2
238	2
240	1
242	2
244	2
246	1
249	1
253	2
255	1
256	1
257	1
258	2
259	2
261	1
263	2
265	1
268	1
269	2
270	1
273	1
277	1
279	1
282	1
283	1
284	1
290	1
292	1
299	1
302	2
303	1
306	1
307	1
310	1
311	1
312	1
314	1
315	1
317	1
318	1
319	2
320	1
321	1
322	1
325	1
328	2
331	1
335	2
338	1
341	1
345	1
350	1
355	1
360	2
370	1
373	1
374	1
390	2
391	1
394	1
397	1
412	1
415	1
426	1
429	1
442	1
444	1
452	1
488	1
494	2
502	1
505	1
509	1
515	1
549	2
551	1
584	1
694	1
701	1
703	1
708	1
710	1
712	1
714	1
809	1
881	1
918	1
962	1
1001	1
1158	1
1264	1
1324	1
1628	1
1848	1
1890	1
1899	1
2227	1
2457	1
2584	1
3405	1
4384	1
9005	1
10491	1
13480	1
23174	1
27760	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
.............................._......................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 6.281GB.  
  PID: 228910;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_counts.txt` (229000)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 6.281GB.  
  PID: 229000;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_plot` (229029)
<pre>
Processing H9_PRO-seq_2
INFO: Found real counts for H9_PRO-seq_2 - Total (M): 19.4822 Unique (M): 18.471612

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.281GB.  
  PID: 229029;	Command: Rscript;	Return code: 0;	Memory used: 0.26GB

> `Library complexity`	QC_hg38/H9_PRO-seq_2_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_2_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 15:09:46) elapsed: 226.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam` (229049)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.281GB.  
  PID: 229049;	Command: samtools;	Return code: 0;	Memory used: 0.017GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv` (229061)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/tmp_H9_PRO-seq_2_PE1_ky4ot7sv'
Processing with 12 cores...
Discarding 100 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270755v1']
Keeping 95 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 6.281GB.  
  PID: 229061;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.01GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	9741100.0	PEPPRO	_RES_

> `PBC2`	9741100.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_unmap.bam` (229296)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.281GB.  
  PID: 229296;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam`

> `Unmapped_reads`	6664276	PEPPRO	_RES_

### Split BAM by strand (06-15 15:10:32) elapsed: 46.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam` (229338)
<pre>
</pre>
Command completed. Elapsed time: 0:01:03. Running peak memory: 6.281GB.  
  PID: 229338;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam` (229395)
<pre>
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 6.281GB.  
  PID: 229395;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 15:12:35) elapsed: 123.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (229457)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 229457;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt` (229458)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.281GB.  
  PID: 229458;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.814GB


> `TSS_coding_score`	33.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt` (229494)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.281GB.  
  PID: 229494;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.821GB


> `TSS_non-coding_score`	12.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt` (229525)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.281GB.  
  PID: 229525;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_2_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_2_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt` (229548,229549,229550,229551)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229548;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 229550;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 229549;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 229551;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt` (229553)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229553;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-15 15:12:56) elapsed: 21.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed` (229555,229556)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.281GB.  
  PID: 229555;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 229556;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed` (229560,229561)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229560;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 229561;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSS_density.bed` (229563,229564,229565,229566)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.281GB.  
  PID: 229563;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 229565;	Command: sort;	Return code: 0;	Memory used: 0.01GB  
  PID: 229564;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 229566;	Command: sort;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_gene_body_density.bed` (229589,229590,229591)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.281GB.  
  PID: 229590;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 229589;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB  
  PID: 229591;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/tmp3wxoh0t5` (229619,229621,229622)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 229619;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 229622;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 229621;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/tmp3wxoh0t5 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0187236) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/tmp3wxoh0t5 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed` (229628)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229628;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	34.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed` (229633)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.281GB.  
  PID: 229633;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Pause index`	QC_hg38/H9_PRO-seq_2_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_2_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed` (229655)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229655;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 15:13:58) elapsed: 62.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam`
18165727.0 7048811

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam`
18165727.0 6572286

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_gene_sort.bed` (229699,229700)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 229699;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 229700;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_gene_coverage.bed` (229703)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.281GB.  
  PID: 229703;	Command: bedtools;	Return code: 0;	Memory used: 0.1GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed.gz` (229731)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229731;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed` (229732)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229732;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 15:15:01) elapsed: 63.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed` (229859)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 229859;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed` (229861,229862,229863,229864)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 229861;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 229862;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 229864;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 229863;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed` (229866)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.281GB.  
  PID: 229866;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_minus_coverage.bed` (229882)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 6.281GB.  
  PID: 229882;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed` (229893,229894,229895,229896)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229893;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 229894;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 229896;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 229895;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed` (229898)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.281GB.  
  PID: 229898;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_minus_coverage.bed` (229913)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.281GB.  
  PID: 229913;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region"` (229925)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229925;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed` (229926,229927,229928,229929)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 229926;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 229928;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 229927;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 229929;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (229932)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.281GB.  
  PID: 229932;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed` (229944)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.281GB.  
  PID: 229944;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR"` (229955)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229955;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed` (229956,229957,229958,229959)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 229956;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 229957;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 229959;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 229958;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed` (229962)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.281GB.  
  PID: 229962;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_minus_coverage.bed` (229973)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 6.281GB.  
  PID: 229973;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR"` (229987)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 229987;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed` (229988,229989,229990,229991)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 229988;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 229989;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 229991;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 229990;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed` (229994)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.281GB.  
  PID: 229994;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_minus_coverage.bed` (230005)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 6.281GB.  
  PID: 230005;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed` (230017,230018,230019,230020)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.281GB.  
  PID: 230017;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 230018;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 230020;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 230019;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed` (230025)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.281GB.  
  PID: 230025;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_minus_coverage.bed` (230038)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.281GB.  
  PID: 230038;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed` (230050,230051,230052,230053)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 230050;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 230052;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 230051;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 230053;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed` (230056)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.281GB.  
  PID: 230056;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_minus_coverage.bed` (230072)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.281GB.  
  PID: 230072;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB


### Plot cFRiF/FRiF (06-15 15:18:00) elapsed: 179.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_2 -z 3099922541 -n 9944138 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed` (230100)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 6.281GB.  
  PID: 230100;	Command: Rscript;	Return code: 0;	Memory used: 0.428GB

> `cFRiF`	QC_hg38/H9_PRO-seq_2_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_2_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_2 -z 3099922541 -n 9944138 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed` (230151)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.281GB.  
  PID: 230151;	Command: Rscript;	Return code: 0;	Memory used: 0.429GB

> `FRiF`	QC_hg38/H9_PRO-seq_2_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_2_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 15:19:05) elapsed: 65.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_exons_sort.bed` (230184,230185)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.281GB.  
  PID: 230184;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 230185;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_introns_sort.bed` (230191,230192,230193)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.281GB.  
  PID: 230191;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 230193;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 230192;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_coverage.bed` (230199)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 6.281GB.  
  PID: 230199;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_coverage.bed` (230226)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 6.281GB.  
  PID: 230226;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.165727)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_rpkm.bed` (230397,230398,230400)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 230397;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 230400;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 230398;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.165727)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_rpkm.bed` (230402,230403,230404)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.281GB.  
  PID: 230402;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 230404;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 230403;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed` (230440,230441,230442)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 230440;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 230442;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 230441;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed --annotate` (230449)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.281GB.  
  PID: 230449;	Command: Rscript;	Return code: 0;	Memory used: 0.317GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_2_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_2_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed` (230469)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.281GB.  
  PID: 230469;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 15:20:09) elapsed: 64.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam` (230477)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.281GB.  
  PID: 230477;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 18165727.0` (230488)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_2_plus_cuttrace_r1on2vue'
Processing with 4 cores...
stdin is empty of data
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:08. Running peak memory: 6.281GB.  
  PID: 230488;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.779GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam` (231768)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.281GB.  
  PID: 231768;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 18165727.0` (231775)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_2_minus_cuttrace_unr_8kvf'
Processing with 4 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr14_KI270723v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270755v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_exact_body_0-mer.bw'
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:08. Running peak memory: 6.281GB.  
  PID: 231775;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.646GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:35:29
*  Total elapsed time (all runs):  4:52:17
*         Peak memory (this run):  6.281 GB
*        Pipeline completed time: 2020-06-15 15:34:39
