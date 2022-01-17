### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_30 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 8 -M 10000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-aj38-16c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/
*  Pipeline started at:   (06-11 17:54:55) elapsed: 1.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.8
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16373 insertions(+), 3522 deletions(-)

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE2.fastq.gz']`
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
*            `recover`:  `False`
*        `sample_name`:  `H9_PRO-seq_30`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE2.fastq.gz

> `File_mb`	732.68	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:54:58) elapsed: 3.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz` (432526)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 432526;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz` (432528)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 432528;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1.fastq` (432529)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 0.002GB.  
  PID: 432529;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2.fastq` (432753)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.002GB.  
  PID: 432753;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	34718962	PEPPRO	_RES_

> `Fastq_reads`	34718962	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz']

### FASTQ processing:  (06-11 17:55:59) elapsed: 61.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30_R1_cutadapt.txt` (432826)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 1.588GB.  
  PID: 432826;	Command: cutadapt;	Return code: 0;	Memory used: 1.588GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq` (432875,432877)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 1.588GB.  
  PID: 432875;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 432877;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	8498087	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq` (432897)
<pre>
Started analysis of H9_PRO-seq_30_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_30_R1_processed.fastq
Analysis complete for H9_PRO-seq_30_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 1.588GB.  
  PID: 432897;	Command: fastqc;	Return code: 0;	Memory used: 0.172GB

> `FastQC report r1`	fastqc/H9_PRO-seq_30_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_noadap.fastq` (432939)
<pre>
[INFO][0m 439732 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 1.588GB.  
  PID: 432939;	Command: seqkit;	Return code: 0;	Memory used: 1.037GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_trimmed.fastq` (433032,433033)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 1.588GB.  
  PID: 433032;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 433033;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	13756701.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	8586429.0	PEPPRO	_RES_

> `Duplicate_reads`	439732.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.4625	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/processed_R1.flag` (433078)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.588GB.  
  PID: 433078;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30_R2_cutadapt.txt` (433080)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 1.588GB.  
  PID: 433080;	Command: cutadapt;	Return code: 0;	Memory used: 1.561GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed.fastq` (433125,433126)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 1.588GB.  
  PID: 433125;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 433126;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	16996174	PEPPRO	_RES_

> `Trim_loss_rate`	51.05	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq` (433150)
<pre>
Started analysis of H9_PRO-seq_30_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_30_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_30_R1_processed.fastq
Analysis complete for H9_PRO-seq_30_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 1.588GB.  
  PID: 433150;	Command: fastqc;	Return code: 0;	Memory used: 0.178GB

> `FastQC report r1`	fastqc/H9_PRO-seq_30_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed.fastq` (433188)
<pre>
Started analysis of H9_PRO-seq_30_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_30_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_30_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_30_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 1.588GB.  
  PID: 433188;	Command: fastqc;	Return code: 0;	Memory used: 0.168GB

> `FastQC report r2`	fastqc/H9_PRO-seq_30_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30.histogram`  

> `fastq_pair -t 31247065 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_noadap.fastq` (433223)
<pre>
Left paired: 8712079		Right paired: 8712079
Left single: 60973		Right single: 804249
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:43. Running peak memory: 1.588GB.  
  PID: 433223;	Command: fastq_pair;	Return code: 0;	Memory used: 1.346GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_noadap.fastq.paired.fq -o H9_PRO-seq_30 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt` (433269)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 1.588GB.  
  PID: 433269;	Command: flash;	Return code: 0;	Memory used: 0.072GB


### Plot adapter insertion distribution (06-11 18:00:21) elapsed: 262.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt -u 8` (433707)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 1.588GB.  
  PID: 433707;	Command: Rscript;	Return code: 0;	Memory used: 0.103GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_30_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_30_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 18:00:32) elapsed: 12.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/cutadapt/H9_PRO-seq_30.hist`

> `Degradation_ratio`	1.0316	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed_dups.fastq` (433744)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 1.588GB.  
  PID: 433744;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/processed_R2.flag` (433750)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.588GB.  
  PID: 433750;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/repaired.flag`  

> `fastq_pair -t 31247065 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed.fastq` (433751)
<pre>
Left paired: 8430812		Right paired: 8430812
Left single: 67275		Right single: 825759
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:52. Running peak memory: 1.67GB.  
  PID: 433751;	Command: fastq_pair;	Return code: 0;	Memory used: 1.67GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq` (433818)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.67GB.  
  PID: 433818;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed.fastq` (433819)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.67GB.  
  PID: 433819;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/repaired.flag` (433820)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.67GB.  
  PID: 433820;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/dups_repaired.flag`  

> `fastq_pair -t 31247065 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed_dups.fastq` (433821)
<pre>
Left paired: 8028930		Right paired: 8028930
Left single: 52170		Right single: 1227641
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 1.873GB.  
  PID: 433821;	Command: fastq_pair;	Return code: 0;	Memory used: 1.873GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_trimmed.fastq` (433884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.873GB.  
  PID: 433884;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed_dups.fastq` (433885)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.873GB.  
  PID: 433885;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/dups_repaired.flag` (433886)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.873GB.  
  PID: 433886;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-11 18:02:32) elapsed: 120.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 18:02:32) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_bt2` (433887)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.873GB.  
  PID: 433887;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq` (433888)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_30 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
8430812 reads; of these:
  8430812 (100.00%) were unpaired; of these:
    7556778 (89.63%) aligned 0 times
    874034 (10.37%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.37% overall alignment rate

> `Aligned_reads_human_rDNA`	1748068.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.29	PEPPRO	_RES_

### Map to human_rDNA (06-11 18:03:37) elapsed: 65.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_dups_bt2` (433981)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.873GB.  
  PID: 433981;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_dups_R2.fq` (433982)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_30 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
874034 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-11 18:04:43) elapsed: 65.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_30 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/tmp6me7f37r -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp.bam` (434060,434061,434062)
<pre>
752657 reads skipped
0 reads lost
7556778 reads; of these:
  7556778 (100.00%) were paired; of these:
    3133266 (41.46%) aligned concordantly 0 times
    3628100 (48.01%) aligned concordantly exactly 1 time
    795412 (10.53%) aligned concordantly >1 times
    ----
    3133266 pairs aligned concordantly 0 times; of these:
      697022 (22.25%) aligned discordantly 1 time
    ----
    2436244 pairs aligned 0 times concordantly or discordantly; of these:
      4872488 mates make up the pairs; of these:
        2000641 (41.06%) aligned 0 times
        1050292 (21.56%) aligned exactly 1 time
        1821555 (37.38%) aligned >1 times
86.76% overall alignment rate
[bam_sort_core] merging from 4 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:21:33. Running peak memory: 3.66GB.  
  PID: 434060;	Command: bowtie2;	Return code: 0;	Memory used: 3.66GB  
  PID: 434061;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 434062;	Command: samtools;	Return code: 0;	Memory used: 0.922GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam` (436644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 3.66GB.  
  PID: 436644;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	13112915	PEPPRO	_RES_

> `QC_filtered_reads`	7664259	PEPPRO	_RES_

> `Aligned_reads`	5448656.5	PEPPRO	_RES_

> `Alignment_rate`	32.06	PEPPRO	_RES_

> `Total_efficiency`	15.69	PEPPRO	_RES_

> `Read_depth`	1.96	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_30 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/tmp6me7f37r -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp_dups.bam` (437428,437434,437435)
<pre>
7276273 reads; of these:
  7276273 (100.00%) were paired; of these:
    2959267 (40.67%) aligned concordantly 0 times
    3546907 (48.75%) aligned concordantly exactly 1 time
    770099 (10.58%) aligned concordantly >1 times
    ----
    2959267 pairs aligned concordantly 0 times; of these:
      681940 (23.04%) aligned discordantly 1 time
    ----
    2277327 pairs aligned 0 times concordantly or discordantly; of these:
      4554654 mates make up the pairs; of these:
        1880372 (41.28%) aligned 0 times
        1025616 (22.52%) aligned exactly 1 time
        1648666 (36.20%) aligned >1 times
87.08% overall alignment rate
[bam_sort_core] merging from 4 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:17:38. Running peak memory: 3.662GB.  
  PID: 437428;	Command: bowtie2;	Return code: 0;	Memory used: 3.662GB  
  PID: 437434;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 437435;	Command: samtools;	Return code: 0;	Memory used: 0.896GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam` (439839)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.662GB.  
  PID: 439839;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (06-11 18:50:44) elapsed: 2761.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq` (439913)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.662GB.  
  PID: 439913;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R1.fq` (439953)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.662GB.  
  PID: 439953;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp.bam` (440011)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.662GB.  
  PID: 440011;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	251868	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam` (440033)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.662GB.  
  PID: 440033;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/chr_sizes.bed` (440047,440048,440049,440050)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 440048;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 440050;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 440047;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 440049;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/chr_sizes.bed -b -@ 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_noMT.bam` (440052)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.662GB.  
  PID: 440052;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam` (440075)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.662GB.  
  PID: 440075;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam` (440076)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.662GB.  
  PID: 440076;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


### Split BAM file (06-11 18:51:57) elapsed: 74.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam` (440088,440089)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 3.662GB.  
  PID: 440088;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 440089;	Command: samtools;	Return code: 0;	Memory used: 1.562GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE2.bam` (440167,440168)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 3.662GB.  
  PID: 440167;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 440168;	Command: samtools;	Return code: 0;	Memory used: 1.253GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp_dups.bam` (440272)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.662GB.  
  PID: 440272;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE1.bam` (440287,440288)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 3.662GB.  
  PID: 440287;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 440288;	Command: samtools;	Return code: 0;	Memory used: 1.562GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE2.bam` (440382,440383)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:15) elapsed: 198.0 _TIME_
Child process 433888 (perl) was already terminated.
Child process 433982 (perl) was already terminated.
[W::bgzf_read_block] EOF marker is absent. The input is probably truncated
Child process 440382 (samtools) terminated after 0 sec.
Child process 440383 (samtools) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/recover.lock.aligned_hg38__H9_PRO-seq_30_dups_PE1.bam
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/recover.lock.aligned_hg38__H9_PRO-seq_30_dups_PE2.bam

### Pipeline failed at:  (06-11 18:55:15) elapsed: 1.0 _TIME_

Total time: 1:00:21
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_30 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 8 -M 10000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba26-28c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/
*  Pipeline started at:   (06-11 19:13:42) elapsed: 4.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.8
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16373 insertions(+), 3522 deletions(-)

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_30`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE2.fastq.gz

> `File_mb`	732.68	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:13:42) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1.fastq`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz']

### FASTQ processing:  (06-11 19:13:42) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/processed_R1.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/processed_R2.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/repaired.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/dups_repaired.flag`  

### Prealignments (06-11 19:13:42) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:13:42) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_bt2` (353510)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 353510;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_bt2

### Map to human_rDNA (06-11 19:13:42) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_dups_bt2` (353511)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 353511;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-11 19:13:42) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam`  

### Compress all unmapped read files (06-11 19:13:42) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R1.fq.gz`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq.gz`  

### Split BAM file (06-11 19:13:43) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE2.bam`  

> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam.bai
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/lock.aligned_hg38__H9_PRO-seq_30_dups_PE1.bam
Overwriting target...
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/lock.aligned_hg38__H9_PRO-seq_30_dups_PE2.bam
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE1.bam` (353516,353517)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 1.555GB.  
  PID: 353516;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 353517;	Command: samtools;	Return code: 0;	Memory used: 1.555GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE2.bam` (353602,353603)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 1.555GB.  
  PID: 353602;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 353603;	Command: samtools;	Return code: 0;	Memory used: 1.246GB


### Calculate library complexity (06-11 19:15:15) elapsed: 92.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE1.bam` (353995)
<pre>
BAM_INPUT
TOTAL READS     = 5818689
COUNTS_SUM      = 5818689
DISTINCT READS  = 5.16937e+06
DISTINCT COUNTS = 156
MAX COUNT       = 10646
COUNTS OF 1     = 4.83223e+06
OBSERVED COUNTS (10647)
1	4832231
2	243921
3	49628
4	18072
5	8620
6	4623
7	2980
8	1916
9	1339
10	992
11	738
12	591
13	472
14	409
15	330
16	256
17	219
18	188
19	161
20	138
21	146
22	112
23	97
24	89
25	71
26	65
27	57
28	57
29	60
30	43
31	55
32	41
33	32
34	35
35	34
36	24
37	27
38	23
39	18
40	26
41	21
42	20
43	18
44	16
45	12
46	12
47	14
48	11
49	14
50	18
51	13
52	12
53	4
54	10
55	9
56	3
57	3
58	8
59	6
60	4
61	6
62	8
63	4
64	6
65	7
66	6
67	6
68	8
69	7
70	6
71	6
72	3
73	4
74	4
76	3
77	1
78	2
79	1
80	3
81	2
82	3
83	4
84	3
85	2
86	4
87	2
88	3
90	1
92	4
93	2
94	1
95	3
96	1
98	2
99	3
101	2
102	1
103	1
104	2
105	2
107	2
108	2
109	1
110	1
111	1
113	2
115	1
116	1
117	1
122	1
125	1
126	1
130	3
132	2
136	1
147	1
149	1
150	2
152	1
154	1
158	1
159	2
160	1
164	1
166	1
169	1
171	1
192	1
214	2
216	1
217	1
221	1
224	1
225	1
274	1
284	1
294	1
343	1
372	1
394	1
423	1
432	1
524	1
611	1
636	1
640	1
767	1
803	1
817	1
1134	1
1436	1
2930	1
3511	1
4334	1
9348	1
10646	1

sample size: 1000000
sample size: 2000000
sample size: 3000000
sample size: 4000000
sample size: 5000000
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 1.555GB.  
  PID: 353995;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE1.bam` (354032)
<pre>
BAM_INPUT
TOTAL READS     = 5818689
DISTINCT READS  = 5.16937e+06
DISTINCT COUNTS = 156
MAX COUNT       = 10646
COUNTS OF 1     = 4.83223e+06
MAX TERMS       = 74
OBSERVED COUNTS (10647)
1	4832231
2	243921
3	49628
4	18072
5	8620
6	4623
7	2980
8	1916
9	1339
10	992
11	738
12	591
13	472
14	409
15	330
16	256
17	219
18	188
19	161
20	138
21	146
22	112
23	97
24	89
25	71
26	65
27	57
28	57
29	60
30	43
31	55
32	41
33	32
34	35
35	34
36	24
37	27
38	23
39	18
40	26
41	21
42	20
43	18
44	16
45	12
46	12
47	14
48	11
49	14
50	18
51	13
52	12
53	4
54	10
55	9
56	3
57	3
58	8
59	6
60	4
61	6
62	8
63	4
64	6
65	7
66	6
67	6
68	8
69	7
70	6
71	6
72	3
73	4
74	4
76	3
77	1
78	2
79	1
80	3
81	2
82	3
83	4
84	3
85	2
86	4
87	2
88	3
90	1
92	4
93	2
94	1
95	3
96	1
98	2
99	3
101	2
102	1
103	1
104	2
105	2
107	2
108	2
109	1
110	1
111	1
113	2
115	1
116	1
117	1
122	1
125	1
126	1
130	3
132	2
136	1
147	1
149	1
150	2
152	1
154	1
158	1
159	2
160	1
164	1
166	1
169	1
171	1
192	1
214	2
216	1
217	1
221	1
224	1
225	1
274	1
284	1
294	1
343	1
372	1
394	1
423	1
432	1
524	1
611	1
636	1
640	1
767	1
803	1
817	1
1134	1
1436	1
2930	1
3511	1
4334	1
9348	1
10646	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
......_...................................._........._..............__...................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 1.555GB.  
  PID: 354032;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_counts.txt` (354174)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 1.555GB.  
  PID: 354174;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_plot` (354187)
<pre>
Processing H9_PRO-seq_30
INFO: Found real counts for H9_PRO-seq_30 - Total (M): 5.844297 Unique (M): 5.818689

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 1.555GB.  
  PID: 354187;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Library complexity`	QC_hg38/H9_PRO-seq_30_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_30_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8535	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 19:16:38) elapsed: 83.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam` (354209)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 1.555GB.  
  PID: 354209;	Command: samtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_bamQC.tsv` (354214)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/tmp_H9_PRO-seq_30_PE1_up4hoyen'
Processing with 8 cores...
Discarding 110 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270725v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 85 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 1.555GB.  
  PID: 354214;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.774GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	2922148.5	PEPPRO	_RES_

> `PBC2`	2922148.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_unmap.bam` (354244)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 1.555GB.  
  PID: 354244;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_temp.bam`

> `Unmapped_reads`	2000641	PEPPRO	_RES_

### Split BAM by strand (06-11 19:16:57) elapsed: 19.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam` (354270)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 1.555GB.  
  PID: 354270;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam` (354305)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 1.555GB.  
  PID: 354305;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 19:17:40) elapsed: 43.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (354324)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.555GB.  
  PID: 354324;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_plus_TssEnrichment.txt` (354325)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 1.555GB.  
  PID: 354325;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.66GB


> `TSS_coding_score`	33.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_minus_TssEnrichment.txt` (354347)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 1.555GB.  
  PID: 354347;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.42GB


> `TSS_non-coding_score`	12.2	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_minus_TssEnrichment.txt` (354369)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 1.555GB.  
  PID: 354369;	Command: Rscript;	Return code: 0;	Memory used: 0.276GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_30_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_30_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt` (354390,354391,354392,354393)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.555GB.  
  PID: 354390;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 354392;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 354391;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 354393;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_keep.txt` (354396)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.555GB.  
  PID: 354396;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-11 19:17:56) elapsed: 16.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_ensembl_tss.bed` (354398,354399)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 1.555GB.  
  PID: 354398;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 354399;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_ensembl_gene_body.bed` (354403,354404)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.555GB.  
  PID: 354403;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 354404;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_TSS_density.bed` (354406,354407,354408,354409)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 1.555GB.  
  PID: 354407;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 354409;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 354406;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 354408;	Command: sort;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_gene_body_density.bed` (354423,354424,354425)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 1.555GB.  
  PID: 354425;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 354423;	Command: bedtools;	Return code: 0;	Memory used: 0.031GB  
  PID: 354424;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/tmppnfs67_d` (354436,354437,354438)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.555GB.  
  PID: 354436;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 354438;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 354437;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/tmppnfs67_d | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:18:19) elapsed: 23.0 _TIME_

Total time: 0:04:41
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/tmppnfs67_d | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_30 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 8 -M 10000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aj40-14c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/
*  Pipeline started at:   (06-14 21:18:39) elapsed: 7.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.8
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16373 insertions(+), 3522 deletions(-)

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_30`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_30pct_PE2.fastq.gz

> `File_mb`	732.68	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:18:40) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R1.fastq`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/H9_PRO-seq_30_R2.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/H9_PRO-seq_30_R2.fastq.gz']

### FASTQ processing:  (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/processed_R1.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/processed_R2.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/repaired.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/fastq/dups_repaired.flag`  

### Prealignments (06-14 21:18:40) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:18:40) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_bt2` (287516)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
  PID: 287516;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_bt2

### Map to human_rDNA (06-14 21:18:40) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/human_rDNA_dups_bt2` (287559)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
  PID: 287559;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_sort.bam`  

### Compress all unmapped read files (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R1.fq.gz`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/prealignments/H9_PRO-seq_30_human_rDNA_unmap_R2.fq.gz`  

### Split BAM file (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE2.bam`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_unmap.bam`  

### Split BAM by strand (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam`  

### Calculate TSS enrichment (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/H9_PRO-seq_30_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_30_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:18:40) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/tmpk8miiqwi` (287629,287640,287642)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 287629;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 287642;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 287640;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/tmpk8miiqwi | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.00674865) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/tmpk8miiqwi > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_pause_index.bed` (287832)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 287832;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	32.88	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_pause_index.bed` (287909)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.204GB.  
  PID: 287909;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Pause index`	QC_hg38/H9_PRO-seq_30_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_30_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_pause_index.bed` (289084)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 289084;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:18:47) elapsed: 7.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam`
5448656.5 2115185

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam`
5448656.5 1971334

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_gene_sort.bed` (289116,289117)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 289116;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 289117;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_gene_coverage.bed` (289119)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 0.204GB.  
  PID: 289119;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/hg38_annotations.bed.gz` (289137)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 289137;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/hg38_annotations.bed` (289138)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 289138;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:19:13) elapsed: 26.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/raw/hg38_annotations.bed` (289147)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 289147;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Enhancer_sort.bed` (289154,289155,289156,289157)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 289154;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 289155;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 289157;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 289156;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Enhancer_plus_coverage.bed` (289159)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 289159;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Enhancer_minus_coverage.bed` (289165)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 289165;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_sort.bed` (289173,289174,289175,289176)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 289173;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 289174;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 289176;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 289175;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_plus_coverage.bed` (289178)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.204GB.  
  PID: 289178;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_minus_coverage.bed` (289184)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 289184;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_Flanking_Region"` (289191)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 289191;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed` (289192,289193,289194,289195)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 289192;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 289194;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 289193;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 289195;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_Flanking_Region_plus_coverage.bed` (289198)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.204GB.  
  PID: 289198;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_Flanking_Region_minus_coverage.bed` (289203)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 289203;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/5_UTR"` (289209)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 289209;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/5_UTR_sort.bed` (289210,289211,289212,289213)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 289210;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 289211;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 289213;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 289212;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_5_UTR_plus_coverage.bed` (289215)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 289215;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_5_UTR_minus_coverage.bed` (289221)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 289221;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/3_UTR"` (289226)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 289226;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/3_UTR_sort.bed` (289227,289228,289229,289230)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 289227;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 289228;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 289230;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB  
  PID: 289229;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_3_UTR_plus_coverage.bed` (289233)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 289233;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_3_UTR_minus_coverage.bed` (289238)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 289238;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Exon_sort.bed` (289469,289470,289471,289472)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 289469;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 289470;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 289472;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 289471;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Exon_plus_coverage.bed` (289476)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.204GB.  
  PID: 289476;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Exon_minus_coverage.bed` (289516)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.204GB.  
  PID: 289516;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Intron_sort.bed` (289527,289528,289529,289530)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 289527;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 289529;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 289528;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 289530;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Intron_plus_coverage.bed` (289536)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.204GB.  
  PID: 289536;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Intron_minus_coverage.bed` (289547)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.204GB.  
  PID: 289547;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB


### Plot cFRiF/FRiF (06-14 21:20:31) elapsed: 78.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_30 -z 3099922541 -n 2983906 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Intron_plus_coverage.bed` (289564)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 0.441GB.  
  PID: 289564;	Command: Rscript;	Return code: 0;	Memory used: 0.441GB

> `cFRiF`	QC_hg38/H9_PRO-seq_30_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_30_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_30 -z 3099922541 -n 2983906 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_Intron_plus_coverage.bed` (289616)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 0.441GB.  
  PID: 289616;	Command: Rscript;	Return code: 0;	Memory used: 0.441GB

> `FRiF`	QC_hg38/H9_PRO-seq_30_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_30_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:21:45) elapsed: 74.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_exons_sort.bed` (289658,289659)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.441GB.  
  PID: 289659;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 289658;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_introns_sort.bed` (289666,289667,289668)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.441GB.  
  PID: 289666;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 289668;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 289667;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exons_coverage.bed` (289674)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 0.441GB.  
  PID: 289674;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_introns_coverage.bed` (292119)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.441GB.  
  PID: 292119;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/5.4486565)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exons_rpkm.bed` (301495,301517,301529)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.441GB.  
  PID: 301495;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 301529;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 301517;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/5.4486565)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_introns_rpkm.bed` (302226,302234,302240)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.441GB.  
  PID: 302226;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 302240;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 302234;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exon_intron_ratios.bed` (303039,303042,303062)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.441GB.  
  PID: 303039;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 303062;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 303042;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exon_intron_ratios.bed --annotate` (303474)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.441GB.  
  PID: 303474;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_30_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_30_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/QC_hg38/H9_PRO-seq_30_exon_intron_ratios.bed` (312794)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.441GB.  
  PID: 312794;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (06-14 21:22:24) elapsed: 40.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam` (312929)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.441GB.  
  PID: 312929;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge --scale 5448656.5` (316552)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_30_plus_cuttrace_fk4hlazk'
Processing with 2 cores...
stdin is empty of data
Discarding 128 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000218v1']
Keeping 67 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrEBV']
Reduce step (merge files)...
Merging 67 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_plus_exact_body_0-mer.bw'
Merging 67 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:50. Running peak memory: 2.218GB.  
  PID: 316552;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.218GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam` (9073)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 2.218GB.  
  PID: 9073;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_minus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge --scale 5448656.5` (9087)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/aligned_hg38/H9_PRO-seq_30_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_30_minus_cuttrace_tb_wms91'
Processing with 2 cores...
stdin is empty of data
Discarding 129 chunk(s) of reads: ['chrM', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1']
Keeping 66 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 66 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_minus_exact_body_0-mer.bw'
Merging 66 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_30/signal_hg38/H9_PRO-seq_30_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:33. Running peak memory: 2.218GB.  
  PID: 9087;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.89GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:23:20
*  Total elapsed time (all runs):  4:11:18
*         Peak memory (this run):  2.2179 GB
*        Pipeline completed time: 2020-06-14 21:41:53
