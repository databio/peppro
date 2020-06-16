### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_60 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_60pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_60pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-33c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/
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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_60pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_60pct_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_60`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_60pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_60pct_PE2.fastq.gz

> `File_mb`	1443.18	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:20) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/H9_PRO-seq_60_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_60pct_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/H9_PRO-seq_60_R1.fastq.gz` (420360)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 420360;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/H9_PRO-seq_60_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/H9_PRO-seq_60_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_60pct_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/H9_PRO-seq_60_R2.fastq.gz` (420363)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 420363;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/H9_PRO-seq_60_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/H9_PRO-seq_60_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1.fastq` (420365)
<pre>
</pre>
Command completed. Elapsed time: 0:02:03. Running peak memory: 0.002GB.  
  PID: 420365;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/H9_PRO-seq_60_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2.fastq` (420489)
<pre>
</pre>
Command completed. Elapsed time: 0:01:31. Running peak memory: 0.002GB.  
  PID: 420489;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	69434420	PEPPRO	_RES_

> `Fastq_reads`	69434420	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/H9_PRO-seq_60_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/H9_PRO-seq_60_R2.fastq.gz']

### FASTQ processing:  (06-15 07:22:05) elapsed: 286.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60_R1_cutadapt.txt` (420986)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 3.008GB.  
  PID: 420986;	Command: cutadapt;	Return code: 0;	Memory used: 3.008GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq` (421131,421132)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 3.008GB.  
  PID: 421131;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 421132;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	16997812	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq` (421189)
<pre>
Started analysis of H9_PRO-seq_60_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_60_R1_processed.fastq
Analysis complete for H9_PRO-seq_60_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.008GB.  
  PID: 421189;	Command: fastqc;	Return code: 0;	Memory used: 0.168GB

> `FastQC report r1`	fastqc/H9_PRO-seq_60_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_noadap.fastq` (421287)
<pre>
[INFO][0m 1324320 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 3.008GB.  
  PID: 421287;	Command: seqkit;	Return code: 0;	Memory used: 1.036GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_trimmed.fastq` (421345,421346)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.008GB.  
  PID: 421345;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 421346;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	27510738.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	17169444.0	PEPPRO	_RES_

> `Duplicate_reads`	1324320.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.4551	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/processed_R1.flag` (421625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.008GB.  
  PID: 421625;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60_R2_cutadapt.txt` (421627)
<pre>
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 3.008GB.  
  PID: 421627;	Command: cutadapt;	Return code: 0;	Memory used: 2.318GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed.fastq` (421709,421710)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.008GB.  
  PID: 421709;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 421710;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	33995624	PEPPRO	_RES_

> `Trim_loss_rate`	51.04	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq` (421745)
<pre>
Started analysis of H9_PRO-seq_60_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_60_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_60_R1_processed.fastq
Analysis complete for H9_PRO-seq_60_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 3.008GB.  
  PID: 421745;	Command: fastqc;	Return code: 0;	Memory used: 0.166GB

> `FastQC report r1`	fastqc/H9_PRO-seq_60_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed.fastq` (421815)
<pre>
Started analysis of H9_PRO-seq_60_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_60_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_60_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_60_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 3.008GB.  
  PID: 421815;	Command: fastqc;	Return code: 0;	Memory used: 0.161GB

> `FastQC report r2`	fastqc/H9_PRO-seq_60_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60.histogram`  

> `fastq_pair -t 62490978 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_noadap.fastq` (421885)
<pre>
Left paired: 17426433		Right paired: 17426433
Left single: 121333		Right single: 1607343
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:57. Running peak memory: 3.94GB.  
  PID: 421885;	Command: fastq_pair;	Return code: 0;	Memory used: 3.94GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_noadap.fastq.paired.fq -o H9_PRO-seq_60 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt` (422539)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 3.94GB.  
  PID: 422539;	Command: flash;	Return code: 0;	Memory used: 0.1GB


### Plot adapter insertion distribution (06-15 07:31:45) elapsed: 580.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt -u 8` (422718)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.94GB.  
  PID: 422718;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_60_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_60_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:31:53) elapsed: 9.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/cutadapt/H9_PRO-seq_60.hist`

> `Degradation_ratio`	1.0317	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed_dups.fastq` (422754)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.94GB.  
  PID: 422754;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/processed_R2.flag` (422789)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 422789;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/repaired.flag`  

> `fastq_pair -t 62490978 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed.fastq` (422790)
<pre>
Left paired: 16863840		Right paired: 16863840
Left single: 133972		Right single: 1650293
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:28. Running peak memory: 3.94GB.  
  PID: 422790;	Command: fastq_pair;	Return code: 0;	Memory used: 3.399GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq` (422883)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 422883;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed.fastq` (422884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 422884;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/repaired.flag` (422886)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 422886;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/dups_repaired.flag`  

> `fastq_pair -t 62490978 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed_dups.fastq` (422888)
<pre>
Left paired: 15645931		Right paired: 15645931
Left single: 100543		Right single: 2868202
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:30. Running peak memory: 3.94GB.  
  PID: 422888;	Command: fastq_pair;	Return code: 0;	Memory used: 3.004GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_trimmed.fastq` (423331)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 423331;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed_dups.fastq` (423334)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 423334;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/dups_repaired.flag` (423335)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 423335;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:35:13) elapsed: 200.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:35:13) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/human_rDNA_bt2` (423337)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 423337;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_R2.fq` (423338)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_60 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
16863840 reads; of these:
  16863840 (100.00%) were unpaired; of these:
    15115450 (89.63%) aligned 0 times
    1748390 (10.37%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.37% overall alignment rate

> `Aligned_reads_human_rDNA`	3496780.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.29	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:38:02) elapsed: 169.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/human_rDNA_dups_bt2` (423572)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 423572;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_dups_R2.fq` (423573)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_60 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/fastq/H9_PRO-seq_60_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
1748390 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:40:22) elapsed: 140.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_60 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/tmpdb81hcn8 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp.bam` (423965,423966,423967)
<pre>
1423235 reads skipped
0 reads lost
15115450 reads; of these:
  15115450 (100.00%) were paired; of these:
    6265719 (41.45%) aligned concordantly 0 times
    7258222 (48.02%) aligned concordantly exactly 1 time
    1591509 (10.53%) aligned concordantly >1 times
    ----
    6265719 pairs aligned concordantly 0 times; of these:
      1393636 (22.24%) aligned discordantly 1 time
    ----
    4872083 pairs aligned 0 times concordantly or discordantly; of these:
      9744166 mates make up the pairs; of these:
        3999373 (41.04%) aligned 0 times
        2100265 (21.55%) aligned exactly 1 time
        3644528 (37.40%) aligned >1 times
86.77% overall alignment rate
[bam_sort_core] merging from 8 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:27:53. Running peak memory: 3.94GB.  
  PID: 423965;	Command: bowtie2;	Return code: 0;	Memory used: 3.768GB  
  PID: 423966;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 423967;	Command: samtools;	Return code: 0;	Memory used: 0.897GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort.bam` (426575)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 3.94GB.  
  PID: 426575;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	26231527	PEPPRO	_RES_

> `QC_filtered_reads`	15331637	PEPPRO	_RES_

> `Aligned_reads`	10899890.0	PEPPRO	_RES_

> `Alignment_rate`	32.06	PEPPRO	_RES_

> `Total_efficiency`	15.7	PEPPRO	_RES_

> `Read_depth`	2.48	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_60 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/tmpdb81hcn8 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp_dups.bam` (428020,428021,428022)
<pre>
14222696 reads; of these:
  14222696 (100.00%) were paired; of these:
    5772579 (40.59%) aligned concordantly 0 times
    6946118 (48.84%) aligned concordantly exactly 1 time
    1503999 (10.57%) aligned concordantly >1 times
    ----
    5772579 pairs aligned concordantly 0 times; of these:
      1336586 (23.15%) aligned discordantly 1 time
    ----
    4435993 pairs aligned 0 times concordantly or discordantly; of these:
      8871986 mates make up the pairs; of these:
        3688993 (41.58%) aligned 0 times
        2009280 (22.65%) aligned exactly 1 time
        3173713 (35.77%) aligned >1 times
87.03% overall alignment rate
[bam_sort_core] merging from 8 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:25:20. Running peak memory: 3.94GB.  
  PID: 428020;	Command: bowtie2;	Return code: 0;	Memory used: 3.755GB  
  PID: 428021;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 428022;	Command: samtools;	Return code: 0;	Memory used: 0.896GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort_dups.bam` (431303)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 3.94GB.  
  PID: 431303;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 08:45:00) elapsed: 3878.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_R1.fq` (431407)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.94GB.  
  PID: 431407;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/prealignments/H9_PRO-seq_60_human_rDNA_unmap_R2.fq` (431647)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.94GB.  
  PID: 431647;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp.bam` (431684)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.94GB.  
  PID: 431684;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	504617	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort.bam` (431743)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.94GB.  
  PID: 431743;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/chr_sizes.bed` (432181,432182,432183,432184)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 432182;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 432184;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 432181;	Command: samtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 432183;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_noMT.bam` (432186)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.94GB.  
  PID: 432186;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort.bam` (432238)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 432238;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort.bam` (432239)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.94GB.  
  PID: 432239;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


### Split BAM file (06-15 08:46:47) elapsed: 108.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam` (432272,432273)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:30. Running peak memory: 3.94GB.  
  PID: 432272;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 432273;	Command: samtools;	Return code: 0;	Memory used: 3.127GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE2.bam` (432399,432400)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:18. Running peak memory: 3.94GB.  
  PID: 432400;	Command: samtools;	Return code: 0;	Memory used: 2.108GB  
  PID: 432399;	Command: samtools;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp_dups.bam` (432882)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.94GB.  
  PID: 432882;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_dups_PE1.bam` (432926,432927)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:29. Running peak memory: 3.94GB.  
  PID: 432926;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 432927;	Command: samtools;	Return code: 0;	Memory used: 3.072GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_dups_PE2.bam` (433155,433156)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:18. Running peak memory: 3.94GB.  
  PID: 433156;	Command: samtools;	Return code: 0;	Memory used: 1.994GB  
  PID: 433155;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


### Calculate library complexity (06-15 08:53:19) elapsed: 392.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_dups_PE1.bam` (433796)
<pre>
BAM_INPUT
TOTAL READS     = 11393512
COUNTS_SUM      = 11393512
DISTINCT READS  = 9.6124e+06
DISTINCT COUNTS = 241
MAX COUNT       = 17651
COUNTS OF 1     = 8.73782e+06
OBSERVED COUNTS (17652)
1	8737820
2	607973
3	135174
4	52518
5	25999
6	14941
7	9257
8	6177
9	4224
10	3156
11	2225
12	1814
13	1399
14	1174
15	964
16	755
17	655
18	539
19	508
20	399
21	362
22	317
23	308
24	251
25	230
26	198
27	167
28	150
29	130
30	155
31	135
32	118
33	135
34	119
35	93
36	93
37	98
38	76
39	66
40	69
41	73
42	54
43	56
44	65
45	34
46	34
47	35
48	37
49	39
50	40
51	28
52	36
53	39
54	29
55	29
56	28
57	32
58	32
59	28
60	21
61	16
62	31
63	26
64	13
65	19
66	17
67	14
68	14
69	8
70	11
71	9
72	13
73	13
74	11
75	10
76	16
77	13
78	7
79	13
80	12
81	4
82	7
83	8
84	11
85	6
86	7
87	10
88	12
89	10
90	7
91	5
92	7
93	2
94	12
95	5
96	5
97	8
98	4
99	4
100	4
101	2
102	1
103	2
104	10
105	7
106	6
107	3
108	4
109	5
110	3
111	3
112	3
113	6
114	2
115	3
116	5
117	2
118	2
119	6
120	3
121	1
123	6
124	3
125	6
126	5
127	2
128	3
129	6
130	1
131	4
132	3
133	2
134	1
135	4
136	1
137	2
139	1
140	6
141	2
142	1
143	4
145	1
147	1
149	1
150	1
153	2
154	1
155	1
156	1
157	3
158	2
159	1
160	1
161	1
162	1
163	1
164	1
165	1
166	1
167	1
168	1
169	1
170	2
172	1
173	1
175	2
177	2
178	1
181	2
182	1
187	1
192	1
193	2
194	1
196	1
197	2
198	1
199	1
200	1
201	1
202	2
203	2
204	1
207	5
209	1
213	1
214	1
215	1
217	1
220	1
222	1
225	1
229	2
231	2
235	1
238	1
250	1
258	1
262	1
263	1
268	1
271	1
288	1
295	1
298	1
305	1
306	1
323	2
325	1
327	1
330	1
341	1
347	1
373	1
423	1
432	2
434	1
436	1
437	2
546	1
551	1
577	1
598	1
637	1
726	1
811	1
815	1
1007	1
1173	1
1196	1
1220	1
1461	1
1489	1
1614	1
2162	1
2790	1
5627	1
6715	1
8463	1
17235	1
17651	1

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
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 3.94GB.  
  PID: 433796;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_dups_PE1.bam` (433884)
<pre>
BAM_INPUT
TOTAL READS     = 11393512
DISTINCT READS  = 9.6124e+06
DISTINCT COUNTS = 241
MAX COUNT       = 17651
COUNTS OF 1     = 8.73782e+06
MAX TERMS       = 100
OBSERVED COUNTS (17652)
1	8737820
2	607973
3	135174
4	52518
5	25999
6	14941
7	9257
8	6177
9	4224
10	3156
11	2225
12	1814
13	1399
14	1174
15	964
16	755
17	655
18	539
19	508
20	399
21	362
22	317
23	308
24	251
25	230
26	198
27	167
28	150
29	130
30	155
31	135
32	118
33	135
34	119
35	93
36	93
37	98
38	76
39	66
40	69
41	73
42	54
43	56
44	65
45	34
46	34
47	35
48	37
49	39
50	40
51	28
52	36
53	39
54	29
55	29
56	28
57	32
58	32
59	28
60	21
61	16
62	31
63	26
64	13
65	19
66	17
67	14
68	14
69	8
70	11
71	9
72	13
73	13
74	11
75	10
76	16
77	13
78	7
79	13
80	12
81	4
82	7
83	8
84	11
85	6
86	7
87	10
88	12
89	10
90	7
91	5
92	7
93	2
94	12
95	5
96	5
97	8
98	4
99	4
100	4
101	2
102	1
103	2
104	10
105	7
106	6
107	3
108	4
109	5
110	3
111	3
112	3
113	6
114	2
115	3
116	5
117	2
118	2
119	6
120	3
121	1
123	6
124	3
125	6
126	5
127	2
128	3
129	6
130	1
131	4
132	3
133	2
134	1
135	4
136	1
137	2
139	1
140	6
141	2
142	1
143	4
145	1
147	1
149	1
150	1
153	2
154	1
155	1
156	1
157	3
158	2
159	1
160	1
161	1
162	1
163	1
164	1
165	1
166	1
167	1
168	1
169	1
170	2
172	1
173	1
175	2
177	2
178	1
181	2
182	1
187	1
192	1
193	2
194	1
196	1
197	2
198	1
199	1
200	1
201	1
202	2
203	2
204	1
207	5
209	1
213	1
214	1
215	1
217	1
220	1
222	1
225	1
229	2
231	2
235	1
238	1
250	1
258	1
262	1
263	1
268	1
271	1
288	1
295	1
298	1
305	1
306	1
323	2
325	1
327	1
330	1
341	1
347	1
373	1
423	1
432	2
434	1
436	1
437	2
546	1
551	1
577	1
598	1
637	1
726	1
811	1
815	1
1007	1
1173	1
1196	1
1220	1
1461	1
1489	1
1614	1
2162	1
2790	1
5627	1
6715	1
8463	1
17235	1
17651	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
....................................................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 3.94GB.  
  PID: 433884;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_counts.txt` (434152)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.94GB.  
  PID: 434152;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_plot` (434189)
<pre>
Processing H9_PRO-seq_60
INFO: Found real counts for H9_PRO-seq_60 - Total (M): 11.690255 Unique (M): 11.393512

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.94GB.  
  PID: 434189;	Command: Rscript;	Return code: 0;	Memory used: 0.318GB

> `Library complexity`	QC_hg38/H9_PRO-seq_60_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_60_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8538	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 08:55:48) elapsed: 149.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam` (434208)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.94GB.  
  PID: 434208;	Command: samtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_bamQC.tsv` (434216)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/tmp_H9_PRO-seq_60_PE1_petsiqdf'
Processing with 12 cores...
Discarding 103 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 92 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.94GB.  
  PID: 434216;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.184GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	5845127.5	PEPPRO	_RES_

> `PBC2`	5845127.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_unmap.bam` (434279)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.94GB.  
  PID: 434279;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_temp.bam`

> `Unmapped_reads`	3999373	PEPPRO	_RES_

### Split BAM by strand (06-15 08:56:18) elapsed: 29.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam` (434333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 3.94GB.  
  PID: 434333;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam` (434368)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 3.94GB.  
  PID: 434368;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 08:57:35) elapsed: 78.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (434425)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434425;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_plus_TssEnrichment.txt` (434427)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.94GB.  
  PID: 434427;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.931GB


> `TSS_coding_score`	33.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_minus_TssEnrichment.txt` (434464)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.94GB.  
  PID: 434464;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.953GB


> `TSS_non-coding_score`	12.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_minus_TssEnrichment.txt` (434495)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.94GB.  
  PID: 434495;	Command: Rscript;	Return code: 0;	Memory used: 0.319GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_60_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_60_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt` (434517,434518,434519,434520)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434517;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 434519;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 434518;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 434520;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_keep.txt` (434522)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434522;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 08:57:52) elapsed: 17.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_ensembl_tss.bed` (434524,434525)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.94GB.  
  PID: 434524;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 434525;	Command: bedtools;	Return code: 0;	Memory used: 0.096GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_ensembl_gene_body.bed` (434528,434529)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434528;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 434529;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_TSS_density.bed` (434531,434532,434533,434534)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.94GB.  
  PID: 434531;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 434533;	Command: sort;	Return code: 0;	Memory used: 0.009GB  
  PID: 434532;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 434534;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_gene_body_density.bed` (434550,434551,434552)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.94GB.  
  PID: 434550;	Command: bedtools;	Return code: 0;	Memory used: 0.071GB  
  PID: 434552;	Command: sort;	Return code: 0;	Memory used: 0.009GB  
  PID: 434551;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/tmp37n822fq` (434570,434571,434572)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434570;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 434572;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 434571;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/tmp37n822fq | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0121273) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/tmp37n822fq > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_pause_index.bed` (434578)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434578;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	33.82	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_pause_index.bed` (434583)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.94GB.  
  PID: 434583;	Command: Rscript;	Return code: 0;	Memory used: 0.319GB

> `Pause index`	QC_hg38/H9_PRO-seq_60_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_60_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_pause_index.bed` (434604)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434604;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 08:58:32) elapsed: 40.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam`
10899890.0 4230357

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam`
10899890.0 3944672

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_gene_sort.bed` (434658,434659)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 434658;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 434659;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_gene_coverage.bed` (434661)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.94GB.  
  PID: 434661;	Command: bedtools;	Return code: 0;	Memory used: 0.071GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/hg38_annotations.bed.gz` (434680)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434680;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/hg38_annotations.bed` (434681)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434681;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 08:59:12) elapsed: 40.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/raw/hg38_annotations.bed` (434689)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 434689;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Enhancer_sort.bed` (434691,434692,434693,434694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 434691;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 434692;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 434694;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB  
  PID: 434693;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Enhancer_plus_coverage.bed` (434697)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.94GB.  
  PID: 434697;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Enhancer_minus_coverage.bed` (434705)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.94GB.  
  PID: 434705;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_sort.bed` (434713,434714,434715,434716)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434713;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 434714;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 434716;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 434715;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_plus_coverage.bed` (434718)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.94GB.  
  PID: 434718;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_minus_coverage.bed` (434727)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.94GB.  
  PID: 434727;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_Flanking_Region"` (434735)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 434735;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed` (434736,434737,434738,434739)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 434736;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 434738;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 434737;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 434739;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_Flanking_Region_plus_coverage.bed` (434742)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.94GB.  
  PID: 434742;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_Flanking_Region_minus_coverage.bed` (434751)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.94GB.  
  PID: 434751;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/5_UTR"` (435030)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 435030;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/5_UTR_sort.bed` (435031,435032,435033,435034)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 435031;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 435032;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 435034;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 435033;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_5_UTR_plus_coverage.bed` (435036)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.94GB.  
  PID: 435036;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_5_UTR_minus_coverage.bed` (435058)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.94GB.  
  PID: 435058;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/3_UTR"` (435091)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 435091;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/3_UTR_sort.bed` (435092,435093,435094,435095)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 435092;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 435093;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 435095;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB  
  PID: 435094;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_3_UTR_plus_coverage.bed` (435098)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.94GB.  
  PID: 435098;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_3_UTR_minus_coverage.bed` (435107)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.94GB.  
  PID: 435107;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Exon_sort.bed` (435131,435132,435133,435134)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.94GB.  
  PID: 435131;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 435132;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 435134;	Command: bedtools;	Return code: 0;	Memory used: 0.17GB  
  PID: 435133;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Exon_plus_coverage.bed` (435138)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.94GB.  
  PID: 435138;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Exon_minus_coverage.bed` (435147)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.94GB.  
  PID: 435147;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Intron_sort.bed` (435155,435156,435157,435158)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 435155;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 435157;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 435156;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 435158;	Command: bedtools;	Return code: 0;	Memory used: 0.082GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Intron_plus_coverage.bed` (435161)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.94GB.  
  PID: 435161;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Intron_minus_coverage.bed` (435183)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.94GB.  
  PID: 435183;	Command: bedtools;	Return code: 0;	Memory used: 0.026GB


### Plot cFRiF/FRiF (06-15 09:01:13) elapsed: 121.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_60 -z 3099922541 -n 5966715 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Intron_plus_coverage.bed` (435206)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.94GB.  
  PID: 435206;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `cFRiF`	QC_hg38/H9_PRO-seq_60_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_60_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_60 -z 3099922541 -n 5966715 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_Intron_plus_coverage.bed` (435273)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.94GB.  
  PID: 435273;	Command: Rscript;	Return code: 0;	Memory used: 0.517GB

> `FRiF`	QC_hg38/H9_PRO-seq_60_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_60_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 09:02:13) elapsed: 60.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_exons_sort.bed` (435305,435306)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.94GB.  
  PID: 435306;	Command: bedtools;	Return code: 0;	Memory used: 0.096GB  
  PID: 435305;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_introns_sort.bed` (435312,435313,435314)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.94GB.  
  PID: 435312;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 435314;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 435313;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exons_coverage.bed` (435320)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.94GB.  
  PID: 435320;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_introns_coverage.bed` (435333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.94GB.  
  PID: 435333;	Command: bedtools;	Return code: 0;	Memory used: 0.061GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/10.89989)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exons_rpkm.bed` (435348,435349,435350)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 435348;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 435350;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 435349;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/10.89989)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_introns_rpkm.bed` (435352,435353,435354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.94GB.  
  PID: 435352;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 435354;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 435353;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exon_intron_ratios.bed` (435358,435359,435360)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 435358;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 435360;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 435359;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exon_intron_ratios.bed --annotate` (435366)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.94GB.  
  PID: 435366;	Command: Rscript;	Return code: 0;	Memory used: 0.32GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_60_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_60_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/QC_hg38/H9_PRO-seq_60_exon_intron_ratios.bed` (435387)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.94GB.  
  PID: 435387;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (06-15 09:02:59) elapsed: 46.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam` (435395)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.94GB.  
  PID: 435395;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 10899890.0` (435400)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_60_plus_cuttrace_rdcn1r3e'
Processing with 4 cores...
stdin is empty of data
Discarding 120 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_GL000009v2_random', 'chr14_KI270725v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 75 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 75 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_plus_exact_body_0-mer.bw'
Merging 75 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:14. Running peak memory: 3.94GB.  
  PID: 435400;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.406GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam` (436902)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.94GB.  
  PID: 436902;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 10899890.0` (436907)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/aligned_hg38/H9_PRO-seq_60_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_60_minus_cuttrace_7wun9_4b'
Processing with 4 cores...
stdin is empty of data
Discarding 116 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270723v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 79 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 79 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_minus_exact_body_0-mer.bw'
Merging 79 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_60/signal_hg38/H9_PRO-seq_60_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:07. Running peak memory: 3.94GB.  
  PID: 436907;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.654GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:58:11
*  Total elapsed time (all runs):  3:23:31
*         Peak memory (this run):  3.9404 GB
*        Pipeline completed time: 2020-06-15 09:15:29
