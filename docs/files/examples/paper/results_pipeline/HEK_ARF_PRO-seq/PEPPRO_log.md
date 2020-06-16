### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name HEK_ARF_PRO-seq --genome hg38 --input /project/shefflab/data//guertin/fastq/SRR8608070_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/SRR8608070_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-32c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/
*  Pipeline started at:   (06-15 07:17:12) elapsed: 1.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/SRR8608070_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/SRR8608070_PE2.fastq.gz']`
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
*        `sample_name`:  `HEK_ARF_PRO-seq`
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

Local input file: /project/shefflab/data//guertin/fastq/SRR8608070_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/SRR8608070_PE2.fastq.gz

> `File_mb`	1302.45	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:13) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/SRR8608070_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R1.fastq.gz` (380686)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 380686;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/SRR8608070_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R2.fastq.gz` (380688)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 380688;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1.fastq` (380695)
<pre>
</pre>
Command completed. Elapsed time: 0:01:44. Running peak memory: 0.002GB.  
  PID: 380695;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2.fastq` (380797)
<pre>
</pre>
Command completed. Elapsed time: 0:00:59. Running peak memory: 0.002GB.  
  PID: 380797;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	53049272	PEPPRO	_RES_

> `Fastq_reads`	53049272	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R2.fastq.gz']

### FASTQ processing:  (06-15 07:20:17) elapsed: 184.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R1_cutadapt.txt` (381168)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 1.701GB.  
  PID: 381168;	Command: cutadapt;	Return code: 0;	Memory used: 1.701GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq` (381268,381269)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 1.701GB.  
  PID: 381268;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 381269;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	16312692	PEPPRO	_RES_

> `Trim_loss_rate`	69.25	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq` (381301)
<pre>
Started analysis of HEK_ARF_PRO-seq_R1_processed.fastq
Approx 5% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 10% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 15% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 20% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 25% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 30% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 35% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 40% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 45% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 50% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 55% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 60% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 65% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 70% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 75% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 80% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 85% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 90% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 95% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Analysis complete for HEK_ARF_PRO-seq_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 1.701GB.  
  PID: 381301;	Command: fastqc;	Return code: 0;	Memory used: 0.157GB

> `FastQC report r1`	fastqc/HEK_ARF_PRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq` (381396)
<pre>
[INFO][0m 778059 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 1.701GB.  
  PID: 381396;	Command: seqkit;	Return code: 0;	Memory used: 1.025GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq` (381446,381448)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 1.701GB.  
  PID: 381446;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 381448;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	21262101.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	8863359.0	PEPPRO	_RES_

> `Duplicate_reads`	778059.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	33.4156	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/processed_R1.flag` (381525)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.701GB.  
  PID: 381525;	Command: touch;	Return code: 0;	Memory used: 0.002GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R2_cutadapt.txt` (381527)
<pre>
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 2.006GB.  
  PID: 381527;	Command: cutadapt;	Return code: 0;	Memory used: 2.006GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq` (381659,381660)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 2.006GB.  
  PID: 381659;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 381660;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	32625384	PEPPRO	_RES_

> `Trim_loss_rate`	38.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq` (381709)
<pre>
Started analysis of HEK_ARF_PRO-seq_R1_processed.fastq
Approx 5% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 10% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 15% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 20% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 25% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 30% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 35% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 40% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 45% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 50% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 55% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 60% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 65% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 70% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 75% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 80% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 85% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 90% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Approx 95% complete for HEK_ARF_PRO-seq_R1_processed.fastq
Analysis complete for HEK_ARF_PRO-seq_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 2.006GB.  
  PID: 381709;	Command: fastqc;	Return code: 0;	Memory used: 0.158GB

> `FastQC report r1`	fastqc/HEK_ARF_PRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq` (381948)
<pre>
Started analysis of HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 5% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 10% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 15% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 20% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 25% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 30% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 35% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 40% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 45% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 50% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 55% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 60% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 65% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 70% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 75% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 80% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 85% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 90% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Approx 95% complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
Analysis complete for HEK_ARF_PRO-seq_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 2.006GB.  
  PID: 381948;	Command: fastqc;	Return code: 0;	Memory used: 0.151GB

> `FastQC report r2`	fastqc/HEK_ARF_PRO-seq_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.histogram`  

> `fastq_pair -t 47744344 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq` (382028)
<pre>
Left paired: 17619354		Right paired: 17619354
Left single: 41924		Right single: 1037248
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:58. Running peak memory: 2.777GB.  
  PID: 382028;	Command: fastq_pair;	Return code: 0;	Memory used: 2.777GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq.paired.fq -o HEK_ARF_PRO-seq -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt` (382093)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 2.777GB.  
  PID: 382093;	Command: flash;	Return code: 0;	Memory used: 0.12GB


### Plot adapter insertion distribution (06-15 07:27:26) elapsed: 429.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt -u 8` (382268)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 2.777GB.  
  PID: 382268;	Command: Rscript;	Return code: 0;	Memory used: 0.21GB

> `Adapter insertion distribution`	cutadapt/HEK_ARF_PRO-seq_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/HEK_ARF_PRO-seq_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	17	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:27:32) elapsed: 6.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`

> `Degradation_ratio`	1.349	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq` (382298)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 2.777GB.  
  PID: 382298;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/processed_R2.flag` (382334)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.777GB.  
  PID: 382334;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/repaired.flag`  

> `fastq_pair -t 47744344 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq` (382335)
<pre>
Left paired: 16265754		Right paired: 16265754
Left single: 46939		Right single: 1128770
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 2.777GB.  
  PID: 382335;	Command: fastq_pair;	Return code: 0;	Memory used: 1.579GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq` (382481)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.777GB.  
  PID: 382481;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq` (382538)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.777GB.  
  PID: 382538;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/repaired.flag` (382539)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.777GB.  
  PID: 382539;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/dups_repaired.flag`  

> `fastq_pair -t 47744344 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq` (382540)
<pre>
Left paired: 15570095		Right paired: 15570095
Left single: 33539		Right single: 1824429
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:05. Running peak memory: 2.777GB.  
  PID: 382540;	Command: fastq_pair;	Return code: 0;	Memory used: 2.327GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq` (382623)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.777GB.  
  PID: 382623;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq` (382624)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.777GB.  
  PID: 382624;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/dups_repaired.flag` (382625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.777GB.  
  PID: 382625;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:29:52) elapsed: 140.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:29:52) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_bt2` (382627)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.777GB.  
  PID: 382627;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq` (382629)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HEK_ARF_PRO-seq -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
16265754 reads; of these:
  16265754 (100.00%) were unpaired; of these:
    14454678 (88.87%) aligned 0 times
    1811076 (11.13%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
11.13% overall alignment rate

> `Aligned_reads_human_rDNA`	3622152.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	11.1	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:31:33) elapsed: 101.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_dups_bt2` (383010)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.777GB.  
  PID: 383010;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_dups_R2.fq` (383011)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HEK_ARF_PRO-seq -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
1811076 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:33:11) elapsed: 98.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id HEK_ARF_PRO-seq -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/tmpyxgj3zkp -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam` (383176,383177,383178)
<pre>
1573721 reads skipped
0 reads lost
14454678 reads; of these:
  14454678 (100.00%) were paired; of these:
    7643921 (52.88%) aligned concordantly 0 times
    5500043 (38.05%) aligned concordantly exactly 1 time
    1310714 (9.07%) aligned concordantly >1 times
    ----
    7643921 pairs aligned concordantly 0 times; of these:
      1654946 (21.65%) aligned discordantly 1 time
    ----
    5988975 pairs aligned 0 times concordantly or discordantly; of these:
      11977950 mates make up the pairs; of these:
        3997391 (33.37%) aligned 0 times
        2623537 (21.90%) aligned exactly 1 time
        5357022 (44.72%) aligned >1 times
86.17% overall alignment rate
[bam_sort_core] merging from 7 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:27:30. Running peak memory: 3.763GB.  
  PID: 383176;	Command: bowtie2;	Return code: 0;	Memory used: 3.763GB  
  PID: 383177;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 383178;	Command: samtools;	Return code: 0;	Memory used: 0.91GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam` (386037)
<pre>
</pre>
Command completed. Elapsed time: 0:00:52. Running peak memory: 3.763GB.  
  PID: 386037;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	24911965	PEPPRO	_RES_

> `QC_filtered_reads`	15309473	PEPPRO	_RES_

> `Aligned_reads`	9602492.0	PEPPRO	_RES_

> `Alignment_rate`	29.43	PEPPRO	_RES_

> `Total_efficiency`	18.1	PEPPRO	_RES_

> `Read_depth`	2.05	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id HEK_ARF_PRO-seq -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/tmpyxgj3zkp -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp_dups.bam` (387068,387074,387075)
<pre>
13996374 reads; of these:
  13996374 (100.00%) were paired; of these:
    7372007 (52.67%) aligned concordantly 0 times
    5351994 (38.24%) aligned concordantly exactly 1 time
    1272373 (9.09%) aligned concordantly >1 times
    ----
    7372007 pairs aligned concordantly 0 times; of these:
      1611786 (21.86%) aligned discordantly 1 time
    ----
    5760221 pairs aligned 0 times concordantly or discordantly; of these:
      11520442 mates make up the pairs; of these:
        3875941 (33.64%) aligned 0 times
        2545638 (22.10%) aligned exactly 1 time
        5098863 (44.26%) aligned >1 times
86.15% overall alignment rate
[bam_sort_core] merging from 7 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:26:14. Running peak memory: 3.763GB.  
  PID: 387068;	Command: bowtie2;	Return code: 0;	Memory used: 3.758GB  
  PID: 387074;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 387075;	Command: samtools;	Return code: 0;	Memory used: 0.908GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam` (390140)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 3.763GB.  
  PID: 390140;	Command: samtools;	Return code: 0;	Memory used: 0.022GB


### Compress all unmapped read files (06-15 08:38:21) elapsed: 3910.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq` (390199)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.763GB.  
  PID: 390199;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R1.fq` (390223)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.763GB.  
  PID: 390223;	Command: pigz;	Return code: 0;	Memory used: 0.013GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam` (390251)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.763GB.  
  PID: 390251;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	252493	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam` (390271)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.763GB.  
  PID: 390271;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/chr_sizes.bed` (390283,390284,390285,390286)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 390284;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 390286;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 390283;	Command: samtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 390285;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_noMT.bam` (390288)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.763GB.  
  PID: 390288;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam` (390321)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 390321;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam` (390323)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.763GB.  
  PID: 390323;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


### Split BAM file (06-15 08:39:51) elapsed: 90.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam` (390337,390338)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 3.763GB.  
  PID: 390338;	Command: samtools;	Return code: 0;	Memory used: 2.324GB  
  PID: 390337;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE2.bam` (390777,390778)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:59. Running peak memory: 3.763GB.  
  PID: 390777;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 390778;	Command: samtools;	Return code: 0;	Memory used: 1.903GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp_dups.bam` (390903)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.763GB.  
  PID: 390903;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE1.bam` (390926,390927)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 3.763GB.  
  PID: 390927;	Command: samtools;	Return code: 0;	Memory used: 2.278GB  
  PID: 390926;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE2.bam` (391029,391031)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:59. Running peak memory: 3.763GB.  
  PID: 391029;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 391031;	Command: samtools;	Return code: 0;	Memory used: 1.902GB


### Calculate library complexity (06-15 08:45:10) elapsed: 319.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE1.bam` (391411)
<pre>
BAM_INPUT
TOTAL READS     = 10681803
COUNTS_SUM      = 10681803
DISTINCT READS  = 9.74978e+06
DISTINCT COUNTS = 191
MAX COUNT       = 10729
COUNTS OF 1     = 9.18474e+06
OBSERVED COUNTS (10730)
1	9184735
2	448142
3	68359
4	20973
5	9240
6	4893
7	3018
8	2018
9	1466
10	1033
11	789
12	698
13	515
14	446
15	342
16	284
17	260
18	236
19	196
20	182
21	148
22	130
23	109
24	104
25	91
26	87
27	79
28	83
29	71
30	56
31	63
32	59
33	34
34	38
35	37
36	30
37	37
38	36
39	29
40	25
41	24
42	30
43	20
44	13
45	16
46	22
47	13
48	24
49	21
50	7
51	10
52	12
53	13
54	9
55	14
56	7
57	12
58	13
59	10
60	12
61	11
62	6
63	5
64	11
65	9
66	5
67	2
68	6
69	5
70	5
71	2
72	6
73	9
74	7
75	5
76	4
77	4
78	3
79	8
80	1
81	2
82	4
83	11
84	2
85	6
86	3
87	3
88	3
89	6
90	1
91	3
93	1
94	4
95	4
96	1
97	1
98	1
99	3
100	2
102	2
103	2
104	2
105	6
106	3
107	1
108	3
109	2
110	3
112	1
113	3
114	4
115	1
117	1
118	1
119	2
120	5
121	1
122	2
124	1
125	2
126	1
127	1
128	3
129	2
133	1
134	1
135	2
138	1
139	3
140	2
142	1
147	1
148	1
153	2
154	1
157	1
159	1
160	1
161	2
162	1
166	1
172	1
173	1
175	1
177	1
182	1
184	1
185	1
186	1
187	1
188	1
192	1
193	1
196	1
198	1
206	1
212	1
216	1
219	1
223	1
233	1
236	1
248	2
267	1
271	1
280	1
286	1
326	1
346	1
351	1
360	1
364	1
396	1
420	1
431	1
461	1
485	1
511	1
536	1
617	1
667	1
771	1
893	1
927	1
974	1
1194	1
1354	1
3429	1
4692	1
7862	1
10729	1

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
</pre>
Command completed. Elapsed time: 0:00:57. Running peak memory: 3.763GB.  
  PID: 391411;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE1.bam` (391464)
<pre>
BAM_INPUT
TOTAL READS     = 10681803
DISTINCT READS  = 9.74978e+06
DISTINCT COUNTS = 191
MAX COUNT       = 10729
COUNTS OF 1     = 9.18474e+06
MAX TERMS       = 90
OBSERVED COUNTS (10730)
1	9184735
2	448142
3	68359
4	20973
5	9240
6	4893
7	3018
8	2018
9	1466
10	1033
11	789
12	698
13	515
14	446
15	342
16	284
17	260
18	236
19	196
20	182
21	148
22	130
23	109
24	104
25	91
26	87
27	79
28	83
29	71
30	56
31	63
32	59
33	34
34	38
35	37
36	30
37	37
38	36
39	29
40	25
41	24
42	30
43	20
44	13
45	16
46	22
47	13
48	24
49	21
50	7
51	10
52	12
53	13
54	9
55	14
56	7
57	12
58	13
59	10
60	12
61	11
62	6
63	5
64	11
65	9
66	5
67	2
68	6
69	5
70	5
71	2
72	6
73	9
74	7
75	5
76	4
77	4
78	3
79	8
80	1
81	2
82	4
83	11
84	2
85	6
86	3
87	3
88	3
89	6
90	1
91	3
93	1
94	4
95	4
96	1
97	1
98	1
99	3
100	2
102	2
103	2
104	2
105	6
106	3
107	1
108	3
109	2
110	3
112	1
113	3
114	4
115	1
117	1
118	1
119	2
120	5
121	1
122	2
124	1
125	2
126	1
127	1
128	3
129	2
133	1
134	1
135	2
138	1
139	3
140	2
142	1
147	1
148	1
153	2
154	1
157	1
159	1
160	1
161	2
162	1
166	1
172	1
173	1
175	1
177	1
182	1
184	1
185	1
186	1
187	1
188	1
192	1
193	1
196	1
198	1
206	1
212	1
216	1
219	1
223	1
233	1
236	1
248	2
267	1
271	1
280	1
286	1
326	1
346	1
351	1
360	1
364	1
396	1
420	1
431	1
461	1
485	1
511	1
536	1
617	1
667	1
771	1
893	1
927	1
974	1
1194	1
1354	1
3429	1
4692	1
7862	1
10729	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
..................................................................................._.................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 3.763GB.  
  PID: 391464;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_counts.txt` (391557)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.763GB.  
  PID: 391557;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_plot` (391575)
<pre>
Processing HEK_ARF_PRO-seq
INFO: Found real counts for HEK_ARF_PRO-seq - Total (M): 10.869678 Unique (M): 10.681803

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.763GB.  
  PID: 391575;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Library complexity`	QC_hg38/HEK_ARF_PRO-seq_preseq_plot.pdf	Library complexity	QC_hg38/HEK_ARF_PRO-seq_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.9164	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 08:47:25) elapsed: 136.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam` (391600)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.763GB.  
  PID: 391600;	Command: samtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_bamQC.tsv` (391608)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/tmp_HEK_ARF_PRO-seq_PE1_hha_ji71'
Processing with 12 cores...
Discarding 96 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270757v1']
Keeping 99 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270508v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.763GB.  
  PID: 391608;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.11GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	5434839.0	PEPPRO	_RES_

> `PBC2`	5434839.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_unmap.bam` (391648)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.763GB.  
  PID: 391648;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam`

> `Unmapped_reads`	3997391	PEPPRO	_RES_

### Split BAM by strand (06-15 08:47:55) elapsed: 30.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam` (391685)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.763GB.  
  PID: 391685;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam` (391720)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.763GB.  
  PID: 391720;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 08:49:03) elapsed: 68.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (391767)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 391767;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_plus_TssEnrichment.txt` (391768)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.763GB.  
  PID: 391768;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.94GB


> `TSS_coding_score`	18.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_minus_TssEnrichment.txt` (391799)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.763GB.  
  PID: 391799;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.989GB


> `TSS_non-coding_score`	2.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_minus_TssEnrichment.txt` (391829)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.763GB.  
  PID: 391829;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/HEK_ARF_PRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/HEK_ARF_PRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt` (391851,391852,391853,391854)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 391851;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 391853;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 391852;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 391854;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt` (391856)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 391856;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 08:49:18) elapsed: 15.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_tss.bed` (391860,391861)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.763GB.  
  PID: 391860;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 391861;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (391865,391866)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 391865;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 391866;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_TSS_density.bed` (391868,391869,391870,391871)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.763GB.  
  PID: 391871;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 391868;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 391870;	Command: sort;	Return code: 0;	Memory used: 0.008GB  
  PID: 391869;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_gene_body_density.bed` (391885,391886,391887)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.763GB.  
  PID: 391885;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB  
  PID: 391887;	Command: sort;	Return code: 0;	Memory used: 0.009GB  
  PID: 391886;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/tmplde2mrc0` (391909,391910,391911)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 391909;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 391911;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 391910;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/tmplde2mrc0 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0100944) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/tmplde2mrc0 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed` (391917)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 391917;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	11.57	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed` (391922)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.763GB.  
  PID: 391922;	Command: Rscript;	Return code: 0;	Memory used: 0.32GB

> `Pause index`	QC_hg38/HEK_ARF_PRO-seq_pause_index.pdf	Pause index	QC_hg38/HEK_ARF_PRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed` (391945)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 391945;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 08:49:56) elapsed: 39.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam`
9602492.0 4000690

> `Plus_FRiP`	0.42	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam`
9602492.0 3788384

> `Minus_FRiP`	0.39	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_gene_sort.bed` (392247,392248)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.763GB.  
  PID: 392247;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 392248;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_gene_coverage.bed` (392251)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.763GB.  
  PID: 392251;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed.gz` (392274)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 392274;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed` (392275)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 392275;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 08:50:35) elapsed: 38.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed` (392283)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.763GB.  
  PID: 392283;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer_sort.bed` (392286,392287,392288,392289)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.763GB.  
  PID: 392286;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 392287;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 392289;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 392288;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_plus_coverage.bed` (392291)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.763GB.  
  PID: 392291;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_minus_coverage.bed` (392299)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.763GB.  
  PID: 392299;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_sort.bed` (392308,392309,392310,392311)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 392308;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 392309;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 392311;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 392310;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_plus_coverage.bed` (392313)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.763GB.  
  PID: 392313;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_minus_coverage.bed` (392321)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.763GB.  
  PID: 392321;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region"` (392328)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 392328;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (392329,392330,392331,392332)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.763GB.  
  PID: 392329;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 392331;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 392330;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 392332;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (392335)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.763GB.  
  PID: 392335;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (392352)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.763GB.  
  PID: 392352;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR"` (392360)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 392360;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR_sort.bed` (392361,392362,392363,392364)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.763GB.  
  PID: 392361;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 392362;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 392364;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 392363;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_plus_coverage.bed` (392367)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.763GB.  
  PID: 392367;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_minus_coverage.bed` (392374)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.763GB.  
  PID: 392374;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR"` (392384)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 392384;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR_sort.bed` (392385,392386,392387,392388)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.763GB.  
  PID: 392385;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 392386;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 392388;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 392387;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_plus_coverage.bed` (392390)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.763GB.  
  PID: 392390;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_minus_coverage.bed` (392398)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.763GB.  
  PID: 392398;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon_sort.bed` (392405,392406,392407,392408)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.763GB.  
  PID: 392405;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 392406;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 392408;	Command: bedtools;	Return code: 0;	Memory used: 0.164GB  
  PID: 392407;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_plus_coverage.bed` (392413)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.763GB.  
  PID: 392413;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_minus_coverage.bed` (392429)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.763GB.  
  PID: 392429;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron_sort.bed` (392437,392438,392439,392440)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.763GB.  
  PID: 392437;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 392439;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 392438;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 392440;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_plus_coverage.bed` (392443)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.763GB.  
  PID: 392443;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_minus_coverage.bed` (392454)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.763GB.  
  PID: 392454;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB


### Plot cFRiF/FRiF (06-15 08:52:25) elapsed: 111.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s HEK_ARF_PRO-seq -z 3099922541 -n 5523870 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_plus_coverage.bed` (392479)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.763GB.  
  PID: 392479;	Command: Rscript;	Return code: 0;	Memory used: 0.45GB

> `cFRiF`	QC_hg38/HEK_ARF_PRO-seq_cFRiF.pdf	cFRiF	QC_hg38/HEK_ARF_PRO-seq_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s HEK_ARF_PRO-seq -z 3099922541 -n 5523870 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_plus_coverage.bed` (392532)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.763GB.  
  PID: 392532;	Command: Rscript;	Return code: 0;	Memory used: 0.467GB

> `FRiF`	QC_hg38/HEK_ARF_PRO-seq_FRiF.pdf	FRiF	QC_hg38/HEK_ARF_PRO-seq_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 08:53:31) elapsed: 65.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_exons_sort.bed` (392576,392577)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.763GB.  
  PID: 392577;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB  
  PID: 392576;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_introns_sort.bed` (392584,392585,392586)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.763GB.  
  PID: 392584;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 392586;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 392585;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_coverage.bed` (392593)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.763GB.  
  PID: 392593;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_coverage.bed` (392608)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.763GB.  
  PID: 392608;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/9.602492)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_rpkm.bed` (392647,392648,392649)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.763GB.  
  PID: 392647;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 392649;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 392648;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/9.602492)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_rpkm.bed` (392651,392652,392653)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.763GB.  
  PID: 392651;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 392653;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 392652;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed` (392656,392657,392658)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 392656;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 392658;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 392657;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed --annotate` (392664)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.763GB.  
  PID: 392664;	Command: Rscript;	Return code: 0;	Memory used: 0.321GB

> `mRNA contamination`	QC_hg38/HEK_ARF_PRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/HEK_ARF_PRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed` (392684)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.763GB.  
  PID: 392684;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 08:54:16) elapsed: 45.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam` (392693)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.763GB.  
  PID: 392693;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 9602492.0` (392697)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_HEK_ARF_PRO-seq_plus_cuttrace_oij6c0gj'
Processing with 4 cores...
Discarding 113 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270752v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 82 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270589v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_exact_body_0-mer.bw'
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:14. Running peak memory: 3.763GB.  
  PID: 392697;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.417GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam` (395252)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.763GB.  
  PID: 395252;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 9602492.0` (395256)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_HEK_ARF_PRO-seq_minus_cuttrace_6_trlrs6'
Processing with 4 cores...
stdin is empty of data
Discarding 107 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_KI270722v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 88 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270508v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 88 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_minus_exact_body_0-mer.bw'
Merging 88 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:07. Running peak memory: 3.763GB.  
  PID: 395256;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.43GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:51:34
*  Total elapsed time (all runs):  3:19:56
*         Peak memory (this run):  3.7633 GB
*        Pipeline completed time: 2020-06-15 09:08:46
