### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_40 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_40pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 8 -M 10000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_40pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aw29-20b
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/
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
*              `cores`:  `8`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_40pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_40pct_PE2.fastq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `10000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `H9_PRO-seq_40`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_40pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_40pct_PE2.fastq.gz

> `File_mb`	971.29	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:19) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/H9_PRO-seq_40_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_40pct_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/H9_PRO-seq_40_R1.fastq.gz` (12548)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 12548;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/H9_PRO-seq_40_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/H9_PRO-seq_40_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_40pct_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/H9_PRO-seq_40_R2.fastq.gz` (12549)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 12549;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/H9_PRO-seq_40_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/H9_PRO-seq_40_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1.fastq` (12550)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 0.002GB.  
  PID: 12550;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/H9_PRO-seq_40_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2.fastq` (12583)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 0.002GB.  
  PID: 12583;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	46292984	PEPPRO	_RES_

> `Fastq_reads`	46292984	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/H9_PRO-seq_40_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/H9_PRO-seq_40_R2.fastq.gz']

### FASTQ processing:  (06-15 07:18:37) elapsed: 78.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40_R1_cutadapt.txt` (12654)
<pre>
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 1.598GB.  
  PID: 12654;	Command: cutadapt;	Return code: 0;	Memory used: 1.598GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq` (12712,12713)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 1.598GB.  
  PID: 12712;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 12713;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	11331720	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq` (12734)
<pre>
Started analysis of H9_PRO-seq_40_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_40_R1_processed.fastq
Analysis complete for H9_PRO-seq_40_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 1.598GB.  
  PID: 12734;	Command: fastqc;	Return code: 0;	Memory used: 0.169GB

> `FastQC report r1`	fastqc/H9_PRO-seq_40_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_noadap.fastq` (13017)
<pre>
[INFO][0m 691509 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 1.598GB.  
  PID: 13017;	Command: seqkit;	Return code: 0;	Memory used: 1.032GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_trimmed.fastq` (13051,13052)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 1.598GB.  
  PID: 13051;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 13052;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	18342898.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	11448553.0	PEPPRO	_RES_

> `Duplicate_reads`	691509.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.4613	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/processed_R1.flag` (13104)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.598GB.  
  PID: 13104;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40_R2_cutadapt.txt` (13106)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 1.598GB.  
  PID: 13106;	Command: cutadapt;	Return code: 0;	Memory used: 1.093GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed.fastq` (13161,13162)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 1.598GB.  
  PID: 13161;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 13162;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	22663440	PEPPRO	_RES_

> `Trim_loss_rate`	51.04	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq` (13218)
<pre>
Started analysis of H9_PRO-seq_40_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_40_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_40_R1_processed.fastq
Analysis complete for H9_PRO-seq_40_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 1.598GB.  
  PID: 13218;	Command: fastqc;	Return code: 0;	Memory used: 0.169GB

> `FastQC report r1`	fastqc/H9_PRO-seq_40_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed.fastq` (13257)
<pre>
Started analysis of H9_PRO-seq_40_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_40_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_40_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_40_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 1.598GB.  
  PID: 13257;	Command: fastqc;	Return code: 0;	Memory used: 0.164GB

> `FastQC report r2`	fastqc/H9_PRO-seq_40_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40.histogram`  

> `fastq_pair -t 41663685 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_noadap.fastq` (13327)
<pre>
Left paired: 11616977		Right paired: 11616977
Left single: 80962		Right single: 1072452
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 2.274GB.  
  PID: 13327;	Command: fastq_pair;	Return code: 0;	Memory used: 2.274GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_noadap.fastq.paired.fq -o H9_PRO-seq_40 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt` (13389)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 2.274GB.  
  PID: 13389;	Command: flash;	Return code: 0;	Memory used: 0.067GB


### Plot adapter insertion distribution (06-15 07:24:55) elapsed: 377.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt -u 8` (13560)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 2.274GB.  
  PID: 13560;	Command: Rscript;	Return code: 0;	Memory used: 0.122GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_40_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_40_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:25:03) elapsed: 8.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/cutadapt/H9_PRO-seq_40.hist`

> `Degradation_ratio`	1.0318	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed_dups.fastq` (13793)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 2.274GB.  
  PID: 13793;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/processed_R2.flag` (13800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.274GB.  
  PID: 13800;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/repaired.flag`  

> `fastq_pair -t 41663685 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed.fastq` (13801)
<pre>
Left paired: 11242330		Right paired: 11242330
Left single: 89390		Right single: 1100903
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 2.361GB.  
  PID: 13801;	Command: fastq_pair;	Return code: 0;	Memory used: 2.361GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq` (13855)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.361GB.  
  PID: 13855;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed.fastq` (13856)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.361GB.  
  PID: 13856;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/repaired.flag` (13857)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.361GB.  
  PID: 13857;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/dups_repaired.flag`  

> `fastq_pair -t 41663685 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed_dups.fastq` (13858)
<pre>
Left paired: 10608843		Right paired: 10608843
Left single: 68392		Right single: 1734390
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 2.428GB.  
  PID: 13858;	Command: fastq_pair;	Return code: 0;	Memory used: 2.428GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_trimmed.fastq` (13906)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.428GB.  
  PID: 13906;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed_dups.fastq` (13907)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.428GB.  
  PID: 13907;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/dups_repaired.flag` (13909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.428GB.  
  PID: 13909;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:27:05) elapsed: 122.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:27:05) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/human_rDNA_bt2` (13910)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.428GB.  
  PID: 13910;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_R2.fq` (13911)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_40 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
11242330 reads; of these:
  11242330 (100.00%) were unpaired; of these:
    10076783 (89.63%) aligned 0 times
    1165547 (10.37%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.37% overall alignment rate

> `Aligned_reads_human_rDNA`	2331094.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.29	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:28:25) elapsed: 81.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/human_rDNA_dups_bt2` (13995)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.428GB.  
  PID: 13995;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_dups_R2.fq` (13996)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_40 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/fastq/H9_PRO-seq_40_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
1165547 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:29:56) elapsed: 90.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_40 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/tmpa7932xtf -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp.bam` (14082,14083,14084)
<pre>
982150 reads skipped
0 reads lost
10076783 reads; of these:
  10076783 (100.00%) were paired; of these:
    4178107 (41.46%) aligned concordantly 0 times
    4837452 (48.01%) aligned concordantly exactly 1 time
    1061224 (10.53%) aligned concordantly >1 times
    ----
    4178107 pairs aligned concordantly 0 times; of these:
      929745 (22.25%) aligned discordantly 1 time
    ----
    3248362 pairs aligned 0 times concordantly or discordantly; of these:
      6496724 mates make up the pairs; of these:
        2667586 (41.06%) aligned 0 times
        1399857 (21.55%) aligned exactly 1 time
        2429281 (37.39%) aligned >1 times
86.76% overall alignment rate
[bam_sort_core] merging from 5 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:25:17. Running peak memory: 3.662GB.  
  PID: 14083;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 14082;	Command: bowtie2;	Return code: 0;	Memory used: 3.662GB  
  PID: 14084;	Command: samtools;	Return code: 0;	Memory used: 0.896GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort.bam` (16987)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 3.662GB.  
  PID: 16987;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	17485980	PEPPRO	_RES_

> `QC_filtered_reads`	10220292	PEPPRO	_RES_

> `Aligned_reads`	7265688.5	PEPPRO	_RES_

> `Alignment_rate`	32.06	PEPPRO	_RES_

> `Total_efficiency`	15.7	PEPPRO	_RES_

> `Read_depth`	2.14	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_40 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/tmpa7932xtf -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp_dups.bam` (17690,17695,17696)
<pre>
9626693 reads; of these:
  9626693 (100.00%) were paired; of these:
    3912255 (40.64%) aligned concordantly 0 times
    4695529 (48.78%) aligned concordantly exactly 1 time
    1018909 (10.58%) aligned concordantly >1 times
    ----
    3912255 pairs aligned concordantly 0 times; of these:
      903759 (23.10%) aligned discordantly 1 time
    ----
    3008496 pairs aligned 0 times concordantly or discordantly; of these:
      6016992 mates make up the pairs; of these:
        2491165 (41.40%) aligned 0 times
        1357608 (22.56%) aligned exactly 1 time
        2168219 (36.03%) aligned >1 times
87.06% overall alignment rate
[bam_sort_core] merging from 5 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:23:59. Running peak memory: 3.662GB.  
  PID: 17690;	Command: bowtie2;	Return code: 0;	Memory used: 3.65GB  
  PID: 17695;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 17696;	Command: samtools;	Return code: 0;	Memory used: 0.896GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort_dups.bam` (20007)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 3.662GB.  
  PID: 20007;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (06-15 08:27:43) elapsed: 3468.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_R1.fq` (20054)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.662GB.  
  PID: 20054;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/prealignments/H9_PRO-seq_40_human_rDNA_unmap_R2.fq` (20078)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.662GB.  
  PID: 20078;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp.bam` (20099)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.662GB.  
  PID: 20099;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	336077	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort.bam` (20117)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.662GB.  
  PID: 20117;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/chr_sizes.bed` (20127,20128,20129,20130)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 20128;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 20130;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 20127;	Command: samtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 20129;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/chr_sizes.bed -b -@ 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_noMT.bam` (20132)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.662GB.  
  PID: 20132;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort.bam` (20158)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 20158;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort.bam` (20159)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.662GB.  
  PID: 20159;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


### Split BAM file (06-15 08:29:06) elapsed: 83.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam` (20169,20170)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 3.662GB.  
  PID: 20169;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20170;	Command: samtools;	Return code: 0;	Memory used: 2.08GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE2.bam` (20511,20512)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 3.662GB.  
  PID: 20511;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20512;	Command: samtools;	Return code: 0;	Memory used: 1.665GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp_dups.bam` (20612)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.662GB.  
  PID: 20612;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_dups_PE1.bam` (20630,20631)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:05. Running peak memory: 3.662GB.  
  PID: 20630;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20631;	Command: samtools;	Return code: 0;	Memory used: 2.058GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_dups_PE2.bam` (20711,20712)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 3.662GB.  
  PID: 20711;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20712;	Command: samtools;	Return code: 0;	Memory used: 1.648GB


### Calculate library complexity (06-15 08:33:47) elapsed: 281.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_dups_PE1.bam` (20786)
<pre>
BAM_INPUT
TOTAL READS     = 7703559
COUNTS_SUM      = 7703559
DISTINCT READS  = 6.71206e+06
DISTINCT COUNTS = 185
MAX COUNT       = 13085
COUNTS OF 1     = 6.20795e+06
OBSERVED COUNTS (13086)
1	6207952
2	358767
3	75714
4	28508
5	14029
6	7631
7	4735
8	3035
9	2128
10	1643
11	1174
12	916
13	690
14	602
15	507
16	401
17	353
18	333
19	270
20	213
21	177
22	163
23	157
24	134
25	146
26	121
27	103
28	88
29	88
30	64
31	55
32	80
33	58
34	53
35	43
36	36
37	46
38	47
39	39
40	41
41	33
42	34
43	27
44	21
45	15
46	24
47	31
48	24
49	22
50	12
51	21
52	16
53	25
54	9
55	14
56	10
57	16
58	16
59	15
60	10
61	7
62	10
63	9
64	10
65	7
66	8
67	4
68	12
69	11
70	8
71	5
72	4
73	4
74	5
75	5
76	3
77	5
78	4
79	3
80	7
81	5
82	4
83	6
84	5
85	3
86	5
87	3
88	10
89	3
90	4
91	7
92	2
93	6
95	2
96	5
97	4
99	2
100	1
101	1
102	2
103	3
104	1
105	4
106	2
107	3
109	1
110	2
111	3
112	1
113	1
114	2
118	1
120	2
121	2
122	1
123	1
124	1
125	3
126	3
128	3
129	2
130	1
132	2
133	2
135	1
137	1
139	2
140	3
141	2
143	1
145	1
148	2
150	1
151	1
152	2
154	1
155	1
158	1
161	1
163	1
165	1
170	2
171	1
174	1
188	1
190	1
201	1
206	1
209	1
211	2
215	1
219	2
221	1
222	1
223	1
231	1
242	1
275	1
283	1
284	1
287	1
288	2
289	1
374	1
375	1
377	1
448	1
473	1
491	1
545	1
572	1
677	1
827	1
831	1
833	1
1005	1
1061	1
1083	1
1493	1
1906	1
3867	1
4608	1
5719	1
12251	1
13085	1

sample size: 1000000
sample size: 2000000
sample size: 3000000
sample size: 4000000
sample size: 5000000
sample size: 6000000
sample size: 7000000
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 3.662GB.  
  PID: 20786;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_dups_PE1.bam` (20822)
<pre>
BAM_INPUT
TOTAL READS     = 7703559
DISTINCT READS  = 6.71206e+06
DISTINCT COUNTS = 185
MAX COUNT       = 13085
COUNTS OF 1     = 6.20795e+06
MAX TERMS       = 92
OBSERVED COUNTS (13086)
1	6207952
2	358767
3	75714
4	28508
5	14029
6	7631
7	4735
8	3035
9	2128
10	1643
11	1174
12	916
13	690
14	602
15	507
16	401
17	353
18	333
19	270
20	213
21	177
22	163
23	157
24	134
25	146
26	121
27	103
28	88
29	88
30	64
31	55
32	80
33	58
34	53
35	43
36	36
37	46
38	47
39	39
40	41
41	33
42	34
43	27
44	21
45	15
46	24
47	31
48	24
49	22
50	12
51	21
52	16
53	25
54	9
55	14
56	10
57	16
58	16
59	15
60	10
61	7
62	10
63	9
64	10
65	7
66	8
67	4
68	12
69	11
70	8
71	5
72	4
73	4
74	5
75	5
76	3
77	5
78	4
79	3
80	7
81	5
82	4
83	6
84	5
85	3
86	5
87	3
88	10
89	3
90	4
91	7
92	2
93	6
95	2
96	5
97	4
99	2
100	1
101	1
102	2
103	3
104	1
105	4
106	2
107	3
109	1
110	2
111	3
112	1
113	1
114	2
118	1
120	2
121	2
122	1
123	1
124	1
125	3
126	3
128	3
129	2
130	1
132	2
133	2
135	1
137	1
139	2
140	3
141	2
143	1
145	1
148	2
150	1
151	1
152	2
154	1
155	1
158	1
161	1
163	1
165	1
170	2
171	1
174	1
188	1
190	1
201	1
206	1
209	1
211	2
215	1
219	2
221	1
222	1
223	1
231	1
242	1
275	1
283	1
284	1
287	1
288	2
289	1
374	1
375	1
377	1
448	1
473	1
491	1
545	1
572	1
677	1
827	1
831	1
833	1
1005	1
1061	1
1083	1
1493	1
1906	1
3867	1
4608	1
5719	1
12251	1
13085	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
..............._.................................._...................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 3.662GB.  
  PID: 20822;	Command: preseq;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_counts.txt` (21058)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.662GB.  
  PID: 21058;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_plot` (21073)
<pre>
Processing H9_PRO-seq_40
INFO: Found real counts for H9_PRO-seq_40 - Total (M): 7.79292 Unique (M): 7.703559

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.662GB.  
  PID: 21073;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Library complexity`	QC_hg38/H9_PRO-seq_40_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_40_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8539	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 08:35:33) elapsed: 106.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam` (21093)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.662GB.  
  PID: 21093;	Command: samtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_bamQC.tsv` (21104)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/tmp_H9_PRO-seq_40_PE1_n6fxgr9j'
Processing with 8 cores...
Discarding 106 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270745v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 89 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.662GB.  
  PID: 21104;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.84GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	3896460.0	PEPPRO	_RES_

> `PBC2`	3896460.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_unmap.bam` (21161)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.662GB.  
  PID: 21161;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_temp.bam`

> `Unmapped_reads`	2667586	PEPPRO	_RES_

### Split BAM by strand (06-15 08:35:58) elapsed: 24.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam` (21190)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 3.662GB.  
  PID: 21190;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam` (21215)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 3.662GB.  
  PID: 21215;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 08:36:55) elapsed: 57.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (21242)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21242;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_plus_TssEnrichment.txt` (21243)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.662GB.  
  PID: 21243;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.567GB


> `TSS_coding_score`	33.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_minus_TssEnrichment.txt` (21266)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.662GB.  
  PID: 21266;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.695GB


> `TSS_non-coding_score`	12.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_minus_TssEnrichment.txt` (21288)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.662GB.  
  PID: 21288;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_40_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_40_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt` (21309,21310,21311,21312)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21309;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 21311;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 21310;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 21312;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_keep.txt` (21314)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21314;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 08:37:12) elapsed: 17.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_ensembl_tss.bed` (21316,21317)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.662GB.  
  PID: 21316;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 21317;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_ensembl_gene_body.bed` (21321,21322)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21321;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 21322;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_TSS_density.bed` (21324,21325,21326,21327)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.662GB.  
  PID: 21326;	Command: sort;	Return code: 0;	Memory used: 0.008GB  
  PID: 21324;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 21327;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 21325;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_gene_body_density.bed` (21339,21340,21341)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.662GB.  
  PID: 21339;	Command: bedtools;	Return code: 0;	Memory used: 0.055GB  
  PID: 21341;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 21340;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/tmp2cnv4nn_` (21356,21357,21358)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21356;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 21358;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 21357;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/tmp2cnv4nn_ | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0086562) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/tmp2cnv4nn_ > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_pause_index.bed` (21364)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21364;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	33.24	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_pause_index.bed` (21369)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.662GB.  
  PID: 21369;	Command: Rscript;	Return code: 0;	Memory used: 0.296GB

> `Pause index`	QC_hg38/H9_PRO-seq_40_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_40_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_pause_index.bed` (21388)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21388;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 08:37:48) elapsed: 36.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam`
7265688.5 2820700

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam`
7265688.5 2629275

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_gene_sort.bed` (21418,21419)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.662GB.  
  PID: 21418;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 21419;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_gene_coverage.bed` (21422)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.662GB.  
  PID: 21422;	Command: bedtools;	Return code: 0;	Memory used: 0.055GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/hg38_annotations.bed.gz` (21437)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21437;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/hg38_annotations.bed` (21438)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21438;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 08:38:19) elapsed: 31.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/raw/hg38_annotations.bed` (21446)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.662GB.  
  PID: 21446;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Enhancer_sort.bed` (21449,21450,21451,21452)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.662GB.  
  PID: 21449;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 21450;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 21452;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 21451;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Enhancer_plus_coverage.bed` (21454)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.662GB.  
  PID: 21454;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Enhancer_minus_coverage.bed` (21461)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.662GB.  
  PID: 21461;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_sort.bed` (21467,21468,21469,21470)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21467;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 21468;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 21470;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 21469;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_plus_coverage.bed` (21472)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.662GB.  
  PID: 21472;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_minus_coverage.bed` (21479)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.662GB.  
  PID: 21479;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_Flanking_Region"` (21485)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21485;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed` (21486,21487,21488,21489)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.662GB.  
  PID: 21486;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 21488;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 21487;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 21489;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_Flanking_Region_plus_coverage.bed` (21492)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.662GB.  
  PID: 21492;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_Flanking_Region_minus_coverage.bed` (21499)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.662GB.  
  PID: 21499;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/5_UTR"` (21505)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21505;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/5_UTR_sort.bed` (21507,21508,21509,21510)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.662GB.  
  PID: 21507;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 21508;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 21510;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 21509;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_5_UTR_plus_coverage.bed` (21512)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.662GB.  
  PID: 21512;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_5_UTR_minus_coverage.bed` (21519)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.662GB.  
  PID: 21519;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/3_UTR"` (21525)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21525;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/3_UTR_sort.bed` (21526,21527,21528,21529)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.662GB.  
  PID: 21526;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 21527;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 21529;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 21528;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_3_UTR_plus_coverage.bed` (21531)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.662GB.  
  PID: 21531;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_3_UTR_minus_coverage.bed` (21539)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.662GB.  
  PID: 21539;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Exon_sort.bed` (21546,21547,21548,21549)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.662GB.  
  PID: 21546;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 21547;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 21549;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 21548;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Exon_plus_coverage.bed` (21554)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.662GB.  
  PID: 21554;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Exon_minus_coverage.bed` (21561)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.662GB.  
  PID: 21561;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Intron_sort.bed` (21568,21569,21570,21571)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.662GB.  
  PID: 21568;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 21570;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 21569;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 21571;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Intron_plus_coverage.bed` (21574)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.662GB.  
  PID: 21574;	Command: bedtools;	Return code: 0;	Memory used: 0.041GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Intron_minus_coverage.bed` (21583)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.662GB.  
  PID: 21583;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB


### Plot cFRiF/FRiF (06-15 08:39:54) elapsed: 95.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_40 -z 3099922541 -n 3978050 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Intron_plus_coverage.bed` (21601)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.662GB.  
  PID: 21601;	Command: Rscript;	Return code: 0;	Memory used: 0.438GB

> `cFRiF`	QC_hg38/H9_PRO-seq_40_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_40_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_40 -z 3099922541 -n 3978050 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_Intron_plus_coverage.bed` (21892)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 3.662GB.  
  PID: 21892;	Command: Rscript;	Return code: 0;	Memory used: 0.438GB

> `FRiF`	QC_hg38/H9_PRO-seq_40_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_40_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 08:41:01) elapsed: 67.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_exons_sort.bed` (21927,21928)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.662GB.  
  PID: 21928;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 21927;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_introns_sort.bed` (21934,21935,21936)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.662GB.  
  PID: 21934;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 21936;	Command: bedtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 21935;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exons_coverage.bed` (21942)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.662GB.  
  PID: 21942;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_introns_coverage.bed` (21953)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.662GB.  
  PID: 21953;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.2656885)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exons_rpkm.bed` (21965,21966,21967)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.662GB.  
  PID: 21965;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 21967;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 21966;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.2656885)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_introns_rpkm.bed` (21970,21971,21972)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.662GB.  
  PID: 21970;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 21972;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 21971;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exon_intron_ratios.bed` (21974,21975,21976)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 21974;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 21976;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 21975;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exon_intron_ratios.bed --annotate` (21983)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.662GB.  
  PID: 21983;	Command: Rscript;	Return code: 0;	Memory used: 0.295GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_40_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_40_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/QC_hg38/H9_PRO-seq_40_exon_intron_ratios.bed` (22003)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 22003;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (06-15 08:41:43) elapsed: 42.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam` (22011)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.662GB.  
  PID: 22011;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge --scale 7265688.5` (22015)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_40_plus_cuttrace_yc274ool'
Processing with 2 cores...
stdin is empty of data
Discarding 124 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 71 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 71 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_plus_exact_body_0-mer.bw'
Merging 71 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:57. Running peak memory: 3.662GB.  
  PID: 22015;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.633GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam` (23614)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.662GB.  
  PID: 23614;	Command: samtools;	Return code: 0;	Memory used: 0.005GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_minus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge --scale 7265688.5` (23618)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/aligned_hg38/H9_PRO-seq_40_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_40_minus_cuttrace_hhgm_zea'
Processing with 2 cores...
stdin is empty of data
Discarding 123 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1']
Keeping 72 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 72 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_minus_exact_body_0-mer.bw'
Merging 72 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_40/signal_hg38/H9_PRO-seq_40_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:51. Running peak memory: 3.662GB.  
  PID: 23618;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.634GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:44:21
*  Total elapsed time (all runs):  3:08:26
*         Peak memory (this run):  3.6624 GB
*        Pipeline completed time: 2020-06-15 09:01:38
