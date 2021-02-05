### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_50 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_50pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_50pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-33c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/
*  Pipeline started at:   (06-15 07:17:19) elapsed: 1.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_50pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_50pct_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_50`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_50pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_50pct_PE2.fastq.gz

> `File_mb`	1207.83	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:20) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/H9_PRO-seq_50_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_50pct_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/H9_PRO-seq_50_R1.fastq.gz` (420359)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 420359;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/H9_PRO-seq_50_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/H9_PRO-seq_50_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_50pct_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/H9_PRO-seq_50_R2.fastq.gz` (420362)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 420362;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/H9_PRO-seq_50_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/H9_PRO-seq_50_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1.fastq` (420366)
<pre>
</pre>
Command completed. Elapsed time: 0:02:04. Running peak memory: 0.002GB.  
  PID: 420366;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/H9_PRO-seq_50_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2.fastq` (420492)
<pre>
</pre>
Command completed. Elapsed time: 0:01:05. Running peak memory: 0.002GB.  
  PID: 420492;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	57859902	PEPPRO	_RES_

> `Fastq_reads`	57859902	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/H9_PRO-seq_50_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/H9_PRO-seq_50_R2.fastq.gz']

### FASTQ processing:  (06-15 07:21:07) elapsed: 228.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50_R1_cutadapt.txt` (420881)
<pre>
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 2.471GB.  
  PID: 420881;	Command: cutadapt;	Return code: 0;	Memory used: 2.471GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq` (420970,420971)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 2.471GB.  
  PID: 420970;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 420971;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	14164231	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq` (421034)
<pre>
Started analysis of H9_PRO-seq_50_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_50_R1_processed.fastq
Analysis complete for H9_PRO-seq_50_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 2.471GB.  
  PID: 421034;	Command: fastqc;	Return code: 0;	Memory used: 0.169GB

> `FastQC report r1`	fastqc/H9_PRO-seq_50_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_noadap.fastq` (421112)
<pre>
[INFO][0m 987746 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 2.471GB.  
  PID: 421112;	Command: seqkit;	Return code: 0;	Memory used: 1.034GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_trimmed.fastq` (421147,421148)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 2.471GB.  
  PID: 421147;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 421148;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	22925851.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	14307341.0	PEPPRO	_RES_

> `Duplicate_reads`	987746.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.4551	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/processed_R1.flag` (421234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.471GB.  
  PID: 421234;	Command: touch;	Return code: 0;	Memory used: 0.002GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50_R2_cutadapt.txt` (421236)
<pre>
</pre>
Command completed. Elapsed time: 0:01:07. Running peak memory: 2.677GB.  
  PID: 421236;	Command: cutadapt;	Return code: 0;	Memory used: 2.677GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed.fastq` (421551,421552)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 2.677GB.  
  PID: 421551;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 421552;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	28328462	PEPPRO	_RES_

> `Trim_loss_rate`	51.04	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq` (421604)
<pre>
Started analysis of H9_PRO-seq_50_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_50_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_50_R1_processed.fastq
Analysis complete for H9_PRO-seq_50_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 2.677GB.  
  PID: 421604;	Command: fastqc;	Return code: 0;	Memory used: 0.167GB

> `FastQC report r1`	fastqc/H9_PRO-seq_50_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed.fastq` (421679)
<pre>
Started analysis of H9_PRO-seq_50_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_50_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_50_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_50_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 2.677GB.  
  PID: 421679;	Command: fastqc;	Return code: 0;	Memory used: 0.164GB

> `FastQC report r2`	fastqc/H9_PRO-seq_50_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50.histogram`  

> `fastq_pair -t 52073911 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_noadap.fastq` (421725)
<pre>
Left paired: 14521420		Right paired: 14521420
Left single: 101190		Right single: 1339977
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:46. Running peak memory: 2.677GB.  
  PID: 421725;	Command: fastq_pair;	Return code: 0;	Memory used: 2.334GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_noadap.fastq.paired.fq -o H9_PRO-seq_50 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt` (421900)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 2.677GB.  
  PID: 421900;	Command: flash;	Return code: 0;	Memory used: 0.116GB


### Plot adapter insertion distribution (06-15 07:28:54) elapsed: 467.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt -u 8` (422100)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 2.677GB.  
  PID: 422100;	Command: Rscript;	Return code: 0;	Memory used: 0.116GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_50_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_50_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:29:04) elapsed: 10.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/cutadapt/H9_PRO-seq_50.hist`

> `Degradation_ratio`	1.0308	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed_dups.fastq` (422134)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 2.677GB.  
  PID: 422134;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/processed_R2.flag` (422144)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.677GB.  
  PID: 422144;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/repaired.flag`  

> `fastq_pair -t 52073911 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed.fastq` (422145)
<pre>
Left paired: 14052561		Right paired: 14052561
Left single: 111670		Right single: 1375764
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 2.677GB.  
  PID: 422145;	Command: fastq_pair;	Return code: 0;	Memory used: 2.174GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq` (422511)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.677GB.  
  PID: 422511;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed.fastq` (422513)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.677GB.  
  PID: 422513;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/repaired.flag` (422514)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.677GB.  
  PID: 422514;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/dups_repaired.flag`  

> `fastq_pair -t 52073911 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed_dups.fastq` (422515)
<pre>
Left paired: 13145891		Right paired: 13145891
Left single: 84515		Right single: 2282434
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:02. Running peak memory: 2.886GB.  
  PID: 422515;	Command: fastq_pair;	Return code: 0;	Memory used: 2.886GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_trimmed.fastq` (422760)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.886GB.  
  PID: 422760;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed_dups.fastq` (422761)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.886GB.  
  PID: 422761;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/dups_repaired.flag` (422762)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.886GB.  
  PID: 422762;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:32:00) elapsed: 176.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:32:00) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/human_rDNA_bt2` (422763)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.886GB.  
  PID: 422763;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_R2.fq` (422764)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_50 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
14052561 reads; of these:
  14052561 (100.00%) were unpaired; of these:
    12594922 (89.63%) aligned 0 times
    1457639 (10.37%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.37% overall alignment rate

> `Aligned_reads_human_rDNA`	2915278.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.29	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:33:34) elapsed: 94.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/human_rDNA_dups_bt2` (422862)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.886GB.  
  PID: 422862;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_dups_R2.fq` (422863)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_50 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/fastq/H9_PRO-seq_50_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
1457639 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:35:04) elapsed: 90.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_50 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/tmpjlrv07a_ -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp.bam` (423321,423322,423323)
<pre>
1205489 reads skipped
0 reads lost
12594922 reads; of these:
  12594922 (100.00%) were paired; of these:
    5220321 (41.45%) aligned concordantly 0 times
    6048178 (48.02%) aligned concordantly exactly 1 time
    1326423 (10.53%) aligned concordantly >1 times
    ----
    5220321 pairs aligned concordantly 0 times; of these:
      1162030 (22.26%) aligned discordantly 1 time
    ----
    4058291 pairs aligned 0 times concordantly or discordantly; of these:
      8116582 mates make up the pairs; of these:
        3332293 (41.06%) aligned 0 times
        1749898 (21.56%) aligned exactly 1 time
        3034391 (37.39%) aligned >1 times
86.77% overall alignment rate
[bam_sort_core] merging from 7 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:22:26. Running peak memory: 3.768GB.  
  PID: 423321;	Command: bowtie2;	Return code: 0;	Memory used: 3.768GB  
  PID: 423322;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 423323;	Command: samtools;	Return code: 0;	Memory used: 0.898GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort.bam` (425491)
<pre>
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 3.768GB.  
  PID: 425491;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	21857551	PEPPRO	_RES_

> `QC_filtered_reads`	12774487	PEPPRO	_RES_

> `Aligned_reads`	9083064.0	PEPPRO	_RES_

> `Alignment_rate`	32.06	PEPPRO	_RES_

> `Total_efficiency`	15.7	PEPPRO	_RES_

> `Read_depth`	2.31	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_50 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/tmpjlrv07a_ -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp_dups.bam` (426431,426436,426437)
<pre>
11940402 reads; of these:
  11940402 (100.00%) were paired; of these:
    4847903 (40.60%) aligned concordantly 0 times
    5829192 (48.82%) aligned concordantly exactly 1 time
    1263307 (10.58%) aligned concordantly >1 times
    ----
    4847903 pairs aligned concordantly 0 times; of these:
      1122097 (23.15%) aligned discordantly 1 time
    ----
    3725806 pairs aligned 0 times concordantly or discordantly; of these:
      7451612 mates make up the pairs; of these:
        3092419 (41.50%) aligned 0 times
        1685392 (22.62%) aligned exactly 1 time
        2673801 (35.88%) aligned >1 times
87.05% overall alignment rate
[bam_sort_core] merging from 6 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:20:02. Running peak memory: 3.768GB.  
  PID: 426431;	Command: bowtie2;	Return code: 0;	Memory used: 3.748GB  
  PID: 426436;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 426437;	Command: samtools;	Return code: 0;	Memory used: 0.9GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort_dups.bam` (428904)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 3.768GB.  
  PID: 428904;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 08:27:18) elapsed: 3135.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_R2.fq` (428957)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.768GB.  
  PID: 428957;	Command: pigz;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/prealignments/H9_PRO-seq_50_human_rDNA_unmap_R1.fq` (428980)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.768GB.  
  PID: 428980;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp.bam` (429003)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.768GB.  
  PID: 429003;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	420662	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort.bam` (429024)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.768GB.  
  PID: 429024;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/chr_sizes.bed` (429037,429038,429039,429040)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 429039;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 429037;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 429040;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 429038;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_noMT.bam` (429042)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.768GB.  
  PID: 429042;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort.bam` (429072)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 429072;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort.bam` (429073)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.768GB.  
  PID: 429073;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


### Split BAM file (06-15 08:28:45) elapsed: 87.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam` (429085,429086)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:17. Running peak memory: 3.768GB.  
  PID: 429086;	Command: samtools;	Return code: 0;	Memory used: 2.163GB  
  PID: 429085;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE2.bam` (429467,429468)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:05. Running peak memory: 3.768GB.  
  PID: 429468;	Command: samtools;	Return code: 0;	Memory used: 2.05GB  
  PID: 429467;	Command: samtools;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp_dups.bam` (429597)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.768GB.  
  PID: 429597;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_dups_PE1.bam` (429618,429619)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:16. Running peak memory: 3.768GB.  
  PID: 429619;	Command: samtools;	Return code: 0;	Memory used: 2.161GB  
  PID: 429618;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_dups_PE2.bam` (429722,429723)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:05. Running peak memory: 3.768GB.  
  PID: 429722;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 429723;	Command: samtools;	Return code: 0;	Memory used: 2.09GB


### Calculate library complexity (06-15 08:34:16) elapsed: 330.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_dups_PE1.bam` (429818)
<pre>
BAM_INPUT
TOTAL READS     = 9562095
COUNTS_SUM      = 9562095
DISTINCT READS  = 8.19046e+06
DISTINCT COUNTS = 211
MAX COUNT       = 15296
COUNTS OF 1     = 7.50584e+06
OBSERVED COUNTS (15297)
1	7505839
2	480908
3	104705
4	39936
5	19742
6	11264
7	6802
8	4520
9	3103
10	2284
11	1772
12	1265
13	1083
14	819
15	713
16	584
17	474
18	384
19	366
20	311
21	289
22	264
23	215
24	189
25	169
26	160
27	146
28	134
29	124
30	105
31	99
32	111
33	93
34	70
35	71
36	61
37	60
38	58
39	60
40	41
41	42
42	38
43	25
44	47
45	39
46	44
47	29
48	43
49	33
50	25
51	33
52	28
53	28
54	18
55	20
56	17
57	21
58	12
59	19
60	8
61	10
62	16
63	14
64	11
65	12
66	18
67	10
68	8
69	16
70	10
71	14
72	11
73	9
74	7
75	8
76	11
77	10
78	4
79	17
80	6
81	6
82	4
83	4
84	5
85	4
86	4
87	2
88	7
89	4
90	4
91	4
92	11
93	2
94	4
95	2
96	3
97	5
98	2
99	7
100	4
101	2
102	4
103	7
104	4
105	6
106	4
107	2
108	5
109	3
110	1
111	4
113	5
114	3
115	5
116	6
117	1
118	2
119	3
120	1
122	2
123	1
126	3
127	3
128	4
129	2
130	1
131	1
133	2
134	2
135	1
137	3
138	2
144	3
145	1
147	1
149	3
151	1
152	2
153	1
156	2
158	1
159	1
160	1
161	1
163	1
164	2
165	2
166	1
167	2
168	1
169	1
170	2
173	1
174	3
177	1
178	1
180	2
181	1
182	2
183	1
188	2
190	1
193	1
198	1
201	1
205	1
208	1
215	1
216	1
224	1
226	1
236	1
243	1
251	1
254	1
259	2
266	1
268	1
272	2
279	1
284	1
289	1
311	1
351	1
355	1
358	1
368	1
369	1
370	2
461	1
471	1
480	1
527	1
544	1
609	1
669	1
689	1
848	1
998	1
1010	1
1025	1
1244	1
1258	1
1353	1
1847	1
2375	1
4751	1
5746	1
7142	1
15010	1
15296	1

sample size: 1000000
sample size: 2000000
sample size: 3000000
sample size: 4000000
sample size: 5000000
sample size: 6000000
sample size: 7000000
sample size: 8000000
sample size: 9000000
</pre>
Command completed. Elapsed time: 0:00:52. Running peak memory: 3.768GB.  
  PID: 429818;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_dups_PE1.bam` (430051)
<pre>
BAM_INPUT
TOTAL READS     = 9562095
DISTINCT READS  = 8.19046e+06
DISTINCT COUNTS = 211
MAX COUNT       = 15296
COUNTS OF 1     = 7.50584e+06
MAX TERMS       = 100
OBSERVED COUNTS (15297)
1	7505839
2	480908
3	104705
4	39936
5	19742
6	11264
7	6802
8	4520
9	3103
10	2284
11	1772
12	1265
13	1083
14	819
15	713
16	584
17	474
18	384
19	366
20	311
21	289
22	264
23	215
24	189
25	169
26	160
27	146
28	134
29	124
30	105
31	99
32	111
33	93
34	70
35	71
36	61
37	60
38	58
39	60
40	41
41	42
42	38
43	25
44	47
45	39
46	44
47	29
48	43
49	33
50	25
51	33
52	28
53	28
54	18
55	20
56	17
57	21
58	12
59	19
60	8
61	10
62	16
63	14
64	11
65	12
66	18
67	10
68	8
69	16
70	10
71	14
72	11
73	9
74	7
75	8
76	11
77	10
78	4
79	17
80	6
81	6
82	4
83	4
84	5
85	4
86	4
87	2
88	7
89	4
90	4
91	4
92	11
93	2
94	4
95	2
96	3
97	5
98	2
99	7
100	4
101	2
102	4
103	7
104	4
105	6
106	4
107	2
108	5
109	3
110	1
111	4
113	5
114	3
115	5
116	6
117	1
118	2
119	3
120	1
122	2
123	1
126	3
127	3
128	4
129	2
130	1
131	1
133	2
134	2
135	1
137	3
138	2
144	3
145	1
147	1
149	3
151	1
152	2
153	1
156	2
158	1
159	1
160	1
161	1
163	1
164	2
165	2
166	1
167	2
168	1
169	1
170	2
173	1
174	3
177	1
178	1
180	2
181	1
182	2
183	1
188	2
190	1
193	1
198	1
201	1
205	1
208	1
215	1
216	1
224	1
226	1
236	1
243	1
251	1
254	1
259	2
266	1
268	1
272	2
279	1
284	1
289	1
311	1
351	1
355	1
358	1
368	1
369	1
370	2
461	1
471	1
480	1
527	1
544	1
609	1
669	1
689	1
848	1
998	1
1010	1
1025	1
1244	1
1258	1
1353	1
1847	1
2375	1
4751	1
5746	1
7142	1
15010	1
15296	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
....................................................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 3.768GB.  
  PID: 430051;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_counts.txt` (430100)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.768GB.  
  PID: 430100;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_plot` (430116)
<pre>
Processing H9_PRO-seq_50
INFO: Found real counts for H9_PRO-seq_50 - Total (M): 9.741802 Unique (M): 9.562095

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.768GB.  
  PID: 430116;	Command: Rscript;	Return code: 0;	Memory used: 0.125GB

> `Library complexity`	QC_hg38/H9_PRO-seq_50_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_50_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8538	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 08:36:24) elapsed: 128.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam` (430140)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.768GB.  
  PID: 430140;	Command: samtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_bamQC.tsv` (430148)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/tmp_H9_PRO-seq_50_PE1_oklojtws'
Processing with 12 cores...
Discarding 105 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270745v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 90 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.768GB.  
  PID: 430148;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.973GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	4870901.0	PEPPRO	_RES_

> `PBC2`	4870901.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_unmap.bam` (430189)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.768GB.  
  PID: 430189;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_temp.bam`

> `Unmapped_reads`	3332293	PEPPRO	_RES_

### Split BAM by strand (06-15 08:36:51) elapsed: 27.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam` (430225)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.768GB.  
  PID: 430225;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam` (430291)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.768GB.  
  PID: 430291;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 08:38:00) elapsed: 69.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (430322)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430322;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_plus_TssEnrichment.txt` (430323)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.768GB.  
  PID: 430323;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.013GB


> `TSS_coding_score`	33.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_minus_TssEnrichment.txt` (430354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.768GB.  
  PID: 430354;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.786GB


> `TSS_non-coding_score`	12.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_minus_TssEnrichment.txt` (430384)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.768GB.  
  PID: 430384;	Command: Rscript;	Return code: 0;	Memory used: 0.32GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_50_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_50_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt` (430405,430407,430408,430409)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430405;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 430408;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 430407;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 430409;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_keep.txt` (430411)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430411;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 08:38:14) elapsed: 14.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_ensembl_tss.bed` (430413,430414)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.768GB.  
  PID: 430413;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 430414;	Command: bedtools;	Return code: 0;	Memory used: 0.096GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_ensembl_gene_body.bed` (430417,430418)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430417;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 430418;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_TSS_density.bed` (430420,430421,430422,430423)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.768GB.  
  PID: 430423;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 430420;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 430422;	Command: sort;	Return code: 0;	Memory used: 0.009GB  
  PID: 430421;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_gene_body_density.bed` (430437,430438,430439)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.768GB.  
  PID: 430437;	Command: bedtools;	Return code: 0;	Memory used: 0.063GB  
  PID: 430439;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 430438;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/tmpakoc64ig` (430455,430456,430457)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430455;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 430457;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 430456;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/tmpakoc64ig | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0104384) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/tmpakoc64ig > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_pause_index.bed` (430463)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430463;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	33.39	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_pause_index.bed` (430468)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.768GB.  
  PID: 430468;	Command: Rscript;	Return code: 0;	Memory used: 0.319GB

> `Pause index`	QC_hg38/H9_PRO-seq_50_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_50_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_pause_index.bed` (430489)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430489;	Command: pigz;	Return code: 0;	Memory used: 0.004GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 08:38:53) elapsed: 39.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam`
9083064.0 3525174

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam`
9083064.0 3287581

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_gene_sort.bed` (430522,430523)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.768GB.  
  PID: 430522;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 430523;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_gene_coverage.bed` (430526)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.768GB.  
  PID: 430526;	Command: bedtools;	Return code: 0;	Memory used: 0.063GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/hg38_annotations.bed.gz` (430541)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430541;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/hg38_annotations.bed` (430542)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430542;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 08:39:29) elapsed: 36.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/raw/hg38_annotations.bed` (430551)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.768GB.  
  PID: 430551;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Enhancer_sort.bed` (430553,430554,430555,430556)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.768GB.  
  PID: 430553;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 430554;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 430556;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 430555;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Enhancer_plus_coverage.bed` (430558)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.768GB.  
  PID: 430558;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Enhancer_minus_coverage.bed` (430566)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.768GB.  
  PID: 430566;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_sort.bed` (430573,430574,430575,430576)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430573;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 430574;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 430576;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 430575;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_plus_coverage.bed` (430578)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.768GB.  
  PID: 430578;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_minus_coverage.bed` (430585)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.768GB.  
  PID: 430585;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_Flanking_Region"` (430594)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430594;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed` (430595,430596,430597,430598)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.768GB.  
  PID: 430595;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 430597;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 430596;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 430598;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_Flanking_Region_plus_coverage.bed` (430601)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.768GB.  
  PID: 430601;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_Flanking_Region_minus_coverage.bed` (430867)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.768GB.  
  PID: 430867;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/5_UTR"` (430878)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430878;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/5_UTR_sort.bed` (430879,430880,430881,430882)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.768GB.  
  PID: 430879;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 430880;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 430882;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 430881;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_5_UTR_plus_coverage.bed` (430884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.768GB.  
  PID: 430884;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_5_UTR_minus_coverage.bed` (430891)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.768GB.  
  PID: 430891;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/3_UTR"` (430900)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 430900;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/3_UTR_sort.bed` (430901,430902,430903,430904)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.768GB.  
  PID: 430901;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 430902;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 430904;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 430903;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_3_UTR_plus_coverage.bed` (430906)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.768GB.  
  PID: 430906;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_3_UTR_minus_coverage.bed` (430914)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.768GB.  
  PID: 430914;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Exon_sort.bed` (430922,430923,430924,430925)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.768GB.  
  PID: 430922;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 430923;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 430925;	Command: bedtools;	Return code: 0;	Memory used: 0.173GB  
  PID: 430924;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Exon_plus_coverage.bed` (430929)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.768GB.  
  PID: 430929;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Exon_minus_coverage.bed` (430937)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.768GB.  
  PID: 430937;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Intron_sort.bed` (430945,430946,430947,430948)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.768GB.  
  PID: 430945;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 430947;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 430946;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 430948;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Intron_plus_coverage.bed` (430951)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.768GB.  
  PID: 430951;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Intron_minus_coverage.bed` (430960)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.768GB.  
  PID: 430960;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB


### Plot cFRiF/FRiF (06-15 08:41:17) elapsed: 108.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_50 -z 3099922541 -n 4972328 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Intron_plus_coverage.bed` (430985)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 3.768GB.  
  PID: 430985;	Command: Rscript;	Return code: 0;	Memory used: 0.475GB

> `cFRiF`	QC_hg38/H9_PRO-seq_50_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_50_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_50 -z 3099922541 -n 4972328 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_Intron_plus_coverage.bed` (431025)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.768GB.  
  PID: 431025;	Command: Rscript;	Return code: 0;	Memory used: 0.475GB

> `FRiF`	QC_hg38/H9_PRO-seq_50_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_50_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 08:42:16) elapsed: 59.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_exons_sort.bed` (431059,431060)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.768GB.  
  PID: 431060;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB  
  PID: 431059;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_introns_sort.bed` (431066,431067,431068)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.768GB.  
  PID: 431066;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 431068;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 431067;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exons_coverage.bed` (431076)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.768GB.  
  PID: 431076;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_introns_coverage.bed` (431088)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.768GB.  
  PID: 431088;	Command: bedtools;	Return code: 0;	Memory used: 0.054GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/9.083064)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exons_rpkm.bed` (431101,431102,431103)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.768GB.  
  PID: 431101;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 431103;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 431102;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/9.083064)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_introns_rpkm.bed` (431105,431106,431107)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.768GB.  
  PID: 431105;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 431107;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 431106;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exon_intron_ratios.bed` (431111,431112,431113)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 431111;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 431113;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 431112;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exon_intron_ratios.bed --annotate` (431119)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.768GB.  
  PID: 431119;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_50_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_50_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/QC_hg38/H9_PRO-seq_50_exon_intron_ratios.bed` (431140)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.768GB.  
  PID: 431140;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 08:43:00) elapsed: 43.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam` (431148)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.768GB.  
  PID: 431148;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 9083064.0` (431152)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_50_plus_cuttrace_q_z2c2xq'
Processing with 4 cores...
stdin is empty of data
Discarding 123 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 72 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 72 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_plus_exact_body_0-mer.bw'
Merging 72 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:15. Running peak memory: 3.768GB.  
  PID: 431152;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.442GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam` (432469)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.768GB.  
  PID: 432469;	Command: samtools;	Return code: 0;	Memory used: 0.005GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 9083064.0` (432473)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/aligned_hg38/H9_PRO-seq_50_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_50_minus_cuttrace_5n418av4'
Processing with 4 cores...
stdin is empty of data
Discarding 119 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270723v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1']
Keeping 76 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 76 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_minus_exact_body_0-mer.bw'
Merging 76 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_50/signal_hg38/H9_PRO-seq_50_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:12. Running peak memory: 3.768GB.  
  PID: 432473;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.438GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:38:18
*  Total elapsed time (all runs):  2:47:09
*         Peak memory (this run):  3.7684 GB
*        Pipeline completed time: 2020-06-15 08:55:36
