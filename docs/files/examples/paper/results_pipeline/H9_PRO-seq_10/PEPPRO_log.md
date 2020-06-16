### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_10 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 4 -M 8000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-ba27-8
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/
*  Pipeline started at:   (06-11 17:53:46) elapsed: 0.0 _TIME_

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
*              `cores`:  `4`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE2.fastq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `8000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `H9_PRO-seq_10`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE2.fastq.gz

> `File_mb`	249.31	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:53:47) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz` (128153)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 128153;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz` (128544)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 128544;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1.fastq` (129047)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.002GB.  
  PID: 129047;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2.fastq` (157719)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.002GB.  
  PID: 157719;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	11560826	PEPPRO	_RES_

> `Fastq_reads`	11560826	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz']

### FASTQ processing:  (06-11 17:54:06) elapsed: 18.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq`  

> `(cutadapt -j 4 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10_R1_cutadapt.txt` (27602)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 0.442GB.  
  PID: 27602;	Command: cutadapt;	Return code: 0;	Memory used: 0.442GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq` (27631,27633)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.442GB.  
  PID: 27631;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 27633;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	2830874	PEPPRO	_RES_

> `Trim_loss_rate`	75.51	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq` (27649)
<pre>
Started analysis of H9_PRO-seq_10_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_10_R1_processed.fastq
Analysis complete for H9_PRO-seq_10_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 0.442GB.  
  PID: 27649;	Command: fastqc;	Return code: 0;	Memory used: 0.159GB

> `FastQC report r1`	fastqc/H9_PRO-seq_10_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_trimmed.fastq`  

> `seqkit rmdup --threads 4 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_noadap.fastq` (27676)
<pre>
[INFO][0m 81668 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.442GB.  
  PID: 27676;	Command: seqkit;	Return code: 0;	Memory used: 0.266GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_trimmed.fastq` (32161,32206)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.442GB.  
  PID: 32161;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 32206;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	4579820.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	2857948.0	PEPPRO	_RES_

> `Duplicate_reads`	81668.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.4419	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/processed_R1.flag` (59403)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.442GB.  
  PID: 59403;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed.fastq`  

> `(cutadapt -j 4 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10_R2_cutadapt.txt` (62872)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 0.442GB.  
  PID: 62872;	Command: cutadapt;	Return code: 0;	Memory used: 0.418GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed.fastq` (172221,172300)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.442GB.  
  PID: 172221;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 172300;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	5661748	PEPPRO	_RES_

> `Trim_loss_rate`	51.03	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq` (189973)
<pre>
Started analysis of H9_PRO-seq_10_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_10_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_10_R1_processed.fastq
Analysis complete for H9_PRO-seq_10_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.442GB.  
  PID: 189973;	Command: fastqc;	Return code: 0;	Memory used: 0.16GB

> `FastQC report r1`	fastqc/H9_PRO-seq_10_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed.fastq` (33844)
<pre>
Started analysis of H9_PRO-seq_10_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_10_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_10_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_10_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 0.442GB.  
  PID: 33844;	Command: fastqc;	Return code: 0;	Memory used: 0.156GB

> `FastQC report r2`	fastqc/H9_PRO-seq_10_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10.histogram`  

> `fastq_pair -t 10404743 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_noadap.fastq` (33873)
<pre>
Left paired: 2902040		Right paired: 2902040
Left single: 20425		Right single: 267392
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 0.455GB.  
  PID: 33873;	Command: fastq_pair;	Return code: 0;	Memory used: 0.455GB


> `flash -q -t 4 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_noadap.fastq.paired.fq -o H9_PRO-seq_10 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt` (33893)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 0.455GB.  
  PID: 33893;	Command: flash;	Return code: 0;	Memory used: 0.07GB


### Plot adapter insertion distribution (06-11 17:56:21) elapsed: 135.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt -u 8` (66723)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.455GB.  
  PID: 66723;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_10_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_10_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 17:56:28) elapsed: 7.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/cutadapt/H9_PRO-seq_10.hist`

> `Degradation_ratio`	1.033	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed_dups.fastq` (127018)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.455GB.  
  PID: 127018;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/processed_R2.flag` (139171)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.455GB.  
  PID: 139171;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/repaired.flag`  

> `fastq_pair -t 10404743 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed.fastq` (140488)
<pre>
Left paired: 2808336		Right paired: 2808336
Left single: 22538		Right single: 274561
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.485GB.  
  PID: 140488;	Command: fastq_pair;	Return code: 0;	Memory used: 0.485GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq` (37490)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.485GB.  
  PID: 37490;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed.fastq` (38555)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.485GB.  
  PID: 38555;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/repaired.flag` (39270)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.485GB.  
  PID: 39270;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/dups_repaired.flag`  

> `fastq_pair -t 10404743 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed_dups.fastq` (39903)
<pre>
Left paired: 2734472		Right paired: 2734472
Left single: 18067		Right single: 348425
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 0.485GB.  
  PID: 39903;	Command: fastq_pair;	Return code: 0;	Memory used: 0.478GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_trimmed.fastq` (136617)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.485GB.  
  PID: 136617;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed_dups.fastq` (137486)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.485GB.  
  PID: 137486;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/dups_repaired.flag` (138385)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.485GB.  
  PID: 138385;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-11 17:56:56) elapsed: 28.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 17:56:56) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_bt2` (138973)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.485GB.  
  PID: 138973;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq` (139468)
<pre>
</pre>

> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_10 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
2808336 reads; of these:
  2808336 (100.00%) were unpaired; of these:
    2516719 (89.62%) aligned 0 times
    291617 (10.38%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.38% overall alignment rate

> `Aligned_reads_human_rDNA`	583234.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.3	PEPPRO	_RES_

### Map to human_rDNA (06-11 17:57:19) elapsed: 23.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_dups_bt2` (146297)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.485GB.  
  PID: 146297;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_dups_R2.fq` (146865)
<pre>
</pre>

> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_10 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
291617 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-11 17:57:44) elapsed: 25.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id H9_PRO-seq_10 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/tmpo1dzl55x -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp.bam` (170089,170091,170092)
<pre>
268755 reads skipped
0 reads lost
2516719 reads; of these:
  2516719 (100.00%) were paired; of these:
    1043013 (41.44%) aligned concordantly 0 times
    1209151 (48.04%) aligned concordantly exactly 1 time
    264555 (10.51%) aligned concordantly >1 times
    ----
    1043013 pairs aligned concordantly 0 times; of these:
      231310 (22.18%) aligned discordantly 1 time
    ----
    811703 pairs aligned 0 times concordantly or discordantly; of these:
      1623406 mates make up the pairs; of these:
        667337 (41.11%) aligned 0 times
        349321 (21.52%) aligned exactly 1 time
        606748 (37.38%) aligned >1 times
86.74% overall alignment rate
[bam_sort_core] merging from 1 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:11:59. Running peak memory: 3.552GB.  
  PID: 170089;	Command: bowtie2;	Return code: 0;	Memory used: 3.552GB  
  PID: 170092;	Command: samtools;	Return code: 0;	Memory used: 0.894GB  
  PID: 170091;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam` (171478)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.552GB.  
  PID: 171478;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	4366101	PEPPRO	_RES_

> `QC_filtered_reads`	2551719	PEPPRO	_RES_

> `Aligned_reads`	1814382.0	PEPPRO	_RES_

> `Alignment_rate`	32.05	PEPPRO	_RES_

> `Total_efficiency`	15.69	PEPPRO	_RES_

> `Read_depth`	1.54	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort_dups.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id H9_PRO-seq_10 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/tmpo1dzl55x -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp_dups.bam` (172092,172093,172094)
<pre>
2465717 reads; of these:
  2465717 (100.00%) were paired; of these:
    1005289 (40.77%) aligned concordantly 0 times
    1199703 (48.66%) aligned concordantly exactly 1 time
    260725 (10.57%) aligned concordantly >1 times
    ----
    1005289 pairs aligned concordantly 0 times; of these:
      229484 (22.83%) aligned discordantly 1 time
    ----
    775805 pairs aligned 0 times concordantly or discordantly; of these:
      1551610 mates make up the pairs; of these:
        636007 (40.99%) aligned 0 times
        345986 (22.30%) aligned exactly 1 time
        569617 (36.71%) aligned >1 times
87.10% overall alignment rate
[bam_sort_core] merging from 1 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:12:01. Running peak memory: 3.552GB.  
  PID: 172092;	Command: bowtie2;	Return code: 0;	Memory used: 3.54GB  
  PID: 172094;	Command: samtools;	Return code: 0;	Memory used: 0.894GB  
  PID: 172093;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort_dups.bam` (173804)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.552GB.  
  PID: 173804;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


### Compress all unmapped read files (06-11 18:24:56) elapsed: 1632.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R1.fq` (173839)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.552GB.  
  PID: 173839;	Command: pigz;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq` (174058)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.552GB.  
  PID: 174058;	Command: pigz;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp.bam` (174076)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.552GB.  
  PID: 174076;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	83974	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam` (174093)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.552GB.  
  PID: 174093;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/chr_sizes.bed` (174100,174101,174102,174103)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.552GB.  
  PID: 174101;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 174103;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 174100;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 174102;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/chr_sizes.bed -b -@ 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_noMT.bam` (174105)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.552GB.  
  PID: 174105;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam` (174125)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.552GB.  
  PID: 174125;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam` (174126)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.552GB.  
  PID: 174126;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Split BAM file (06-11 18:25:26) elapsed: 30.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam | samtools sort - -@ 4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam` (174133,174134)
<pre>
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.552GB.  
  PID: 174133;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 174134;	Command: samtools;	Return code: 0;	Memory used: 0.545GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam | samtools sort - -@ 4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE2.bam` (174188,174189)
<pre>
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.552GB.  
  PID: 174188;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 174189;	Command: samtools;	Return code: 0;	Memory used: 0.433GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp_dups.bam` (174253)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.552GB.  
  PID: 174253;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort_dups.bam | samtools sort - -@ 4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_dups_PE1.bam` (174266,174267)
<pre>
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.552GB.  
  PID: 174266;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 174267;	Command: samtools;	Return code: 0;	Memory used: 0.535GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort_dups.bam | samtools sort - -@ 4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_dups_PE2.bam` (174317,174318)
<pre>
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.552GB.  
  PID: 174317;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 174318;	Command: samtools;	Return code: 0;	Memory used: 0.461GB


### Calculate library complexity (06-11 18:26:46) elapsed: 80.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_dups_PE1.bam` (174366)
<pre>
BAM_INPUT
TOTAL READS     = 1966383
COUNTS_SUM      = 1966383
DISTINCT READS  = 1.84264e+06
DISTINCT COUNTS = 84
MAX COUNT       = 4335
COUNTS OF 1     = 1.77576e+06
OBSERVED COUNTS (4336)
1	1775755
2	50849
3	8958
4	2905
5	1374
6	805
7	492
8	321
9	242
10	139
11	154
12	84
13	81
14	60
15	57
16	38
17	37
18	31
19	28
20	25
21	17
22	26
23	13
24	9
25	9
26	8
27	7
28	10
29	7
30	8
31	10
32	5
33	2
34	4
35	2
36	1
37	4
38	3
39	5
40	1
41	5
42	4
43	2
45	1
46	1
49	1
50	2
51	1
54	2
55	1
56	2
58	1
59	3
64	1
65	2
67	1
69	1
74	1
75	1
77	1
78	1
81	1
82	1
89	1
95	1
116	1
123	1
153	1
159	1
161	1
170	1
206	1
218	1
245	1
277	1
286	1
289	1
401	1
513	1
991	1
1224	1
1487	1
3220	1
4335	1

sample size: 1000000
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 3.552GB.  
  PID: 174366;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_dups_PE1.bam` (174388)
<pre>
BAM_INPUT
TOTAL READS     = 1966383
DISTINCT READS  = 1.84264e+06
DISTINCT COUNTS = 84
MAX COUNT       = 4335
COUNTS OF 1     = 1.77576e+06
MAX TERMS       = 42
OBSERVED COUNTS (4336)
1	1775755
2	50849
3	8958
4	2905
5	1374
6	805
7	492
8	321
9	242
10	139
11	154
12	84
13	81
14	60
15	57
16	38
17	37
18	31
19	28
20	25
21	17
22	26
23	13
24	9
25	9
26	8
27	7
28	10
29	7
30	8
31	10
32	5
33	2
34	4
35	2
36	1
37	4
38	3
39	5
40	1
41	5
42	4
43	2
45	1
46	1
49	1
50	2
51	1
54	2
55	1
56	2
58	1
59	3
64	1
65	2
67	1
69	1
74	1
75	1
77	1
78	1
81	1
82	1
89	1
95	1
116	1
123	1
153	1
159	1
161	1
170	1
206	1
218	1
245	1
277	1
286	1
289	1
401	1
513	1
991	1
1224	1
1487	1
3220	1
4335	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
.................................................._..............................._...................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.552GB.  
  PID: 174388;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_counts.txt` (174413)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.552GB.  
  PID: 174413;	Command: echo;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_plot` (174425)
<pre>
Processing H9_PRO-seq_10
INFO: Found real counts for H9_PRO-seq_10 - Total (M): 1.945718 Unique (M): 1.966383

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.552GB.  
  PID: 174425;	Command: Rscript;	Return code: 0;	Memory used: 0.185GB

> `Library complexity`	QC_hg38/H9_PRO-seq_10_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_10_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8545	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 18:27:22) elapsed: 36.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam` (174454)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.552GB.  
  PID: 174454;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam -c 4 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_bamQC.tsv` (174458)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/tmp_H9_PRO-seq_10_PE1_o9vn6o75'
Processing with 4 cores...
Discarding 128 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000224v1', 'chrUn_GL000226v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 67 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270730v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.552GB.  
  PID: 174458;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.36GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	972859.0	PEPPRO	_RES_

> `PBC2`	972859.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_unmap.bam`  

> `samtools view -b -@ 4 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_unmap.bam` (174485)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.552GB.  
  PID: 174485;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_temp.bam`

> `Unmapped_reads`	667337	PEPPRO	_RES_

### Split BAM by strand (06-11 18:27:34) elapsed: 12.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam` (174507)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.552GB.  
  PID: 174507;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam` (174522)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.552GB.  
  PID: 174522;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 18:27:49) elapsed: 15.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (174539)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.552GB.  
  PID: 174539;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_plus_TssEnrichment.txt` (174541)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.552GB.  
  PID: 174541;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.224GB


> `TSS_coding_score`	33.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_minus_TssEnrichment.txt` (174558)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.552GB.  
  PID: 174558;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.227GB


> `TSS_non-coding_score`	12.2	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_minus_TssEnrichment.txt` (174580)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.552GB.  
  PID: 174580;	Command: Rscript;	Return code: 0;	Memory used: 0.275GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_10_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_10_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt` (174619,174620,174621,174622)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.552GB.  
  PID: 174619;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 174621;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 174620;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 174622;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_keep.txt` (174624)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.552GB.  
  PID: 174624;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-11 18:28:04) elapsed: 15.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_ensembl_tss.bed` (174626,174627)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.552GB.  
  PID: 174626;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 174627;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_ensembl_gene_body.bed` (174635,174636)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.552GB.  
  PID: 174635;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 174636;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_TSS_density.bed` (174638,174639,174640,174641)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.552GB.  
  PID: 174638;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 174640;	Command: sort;	Return code: 0;	Memory used: 0.007GB  
  PID: 174639;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 174641;	Command: sort;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_gene_body_density.bed` (174652,174653,174654)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.552GB.  
  PID: 174652;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 174654;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 174653;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_pause_index.bed`  
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:             ^ syntax error
awk: cmd. line:1: { if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}
awk: cmd. line:1:                                                                                                                      ^ syntax error

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == ){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/tmpqe2mrkhv` (174664,174665,174666)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:15) elapsed: 1631.0 _TIME_
Child process 139468 (perl) was already terminated.
Child process 146865 (perl) was already terminated.
Child process 174664 (join) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/recover.lock.QC_hg38__H9_PRO-seq_10_pause_index.bed

### Pipeline failed at:  (06-11 18:55:15) elapsed: 0.0 _TIME_

Total time: 1:01:29
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_10 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 4 -M 8000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-32c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/
*  Pipeline started at:   (06-11 19:13:06) elapsed: 0.0 _TIME_

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
*              `cores`:  `4`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE2.fastq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `8000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `H9_PRO-seq_10`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE2.fastq.gz

> `File_mb`	249.31	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:13:07) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1.fastq`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz']

### FASTQ processing:  (06-11 19:13:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/processed_R1.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/processed_R2.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/repaired.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/dups_repaired.flag`  

### Prealignments (06-11 19:13:07) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:13:07) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_bt2` (83103)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 83103;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_bt2

### Map to human_rDNA (06-11 19:13:07) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_dups_bt2` (83104)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 83104;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-11 19:13:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam`  

### Compress all unmapped read files (06-11 19:13:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R1.fq.gz`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq.gz`  

### Split BAM file (06-11 19:13:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE2.bam`  

### Calculate NRF, PBC1, and PBC2 (06-11 19:13:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_unmap.bam`  

### Split BAM by strand (06-11 19:13:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam`  

### Calculate TSS enrichment (06-11 19:13:07) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/H9_PRO-seq_10_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_10_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-11 19:13:07) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_gene_body_density.bed`  
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/lock.QC_hg38__H9_PRO-seq_10_pause_index.bed
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/tmp5czd4m_l` (83111,83112,83114)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 83111;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 83114;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 83112;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/tmp5czd4m_l | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:13:08) elapsed: 0.0 _TIME_

Total time: 0:00:02
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/tmp5czd4m_l | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_10 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 4 -M 8000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aj37-17c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/
*  Pipeline started at:   (06-14 21:18:39) elapsed: 5.0 _TIME_

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
*              `cores`:  `4`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE2.fastq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `8000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `H9_PRO-seq_10`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_10pct_PE2.fastq.gz

> `File_mb`	249.31	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:18:40) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R1.fastq`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/H9_PRO-seq_10_R2.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/H9_PRO-seq_10_R2.fastq.gz']

### FASTQ processing:  (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/processed_R1.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/processed_R2.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/repaired.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/fastq/dups_repaired.flag`  

### Prealignments (06-14 21:18:40) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:18:40) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_bt2` (315116)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 315116;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_bt2

### Map to human_rDNA (06-14 21:18:40) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/human_rDNA_dups_bt2` (315213)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 315213;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_sort.bam`  

### Compress all unmapped read files (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R1.fq.gz`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/prealignments/H9_PRO-seq_10_human_rDNA_unmap_R2.fq.gz`  

### Split BAM file (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE2.bam`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_unmap.bam`  

### Split BAM by strand (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam`  

### Calculate TSS enrichment (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/H9_PRO-seq_10_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_10_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:18:40) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/tmplttv258v` (315424,315431,315443)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 315424;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 315443;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 315431;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/tmplttv258v | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.00261755) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/tmplttv258v > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_pause_index.bed` (316092)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 316092;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	33.58	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_pause_index.bed` (316301)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.204GB.  
  PID: 316301;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Pause index`	QC_hg38/H9_PRO-seq_10_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_10_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_pause_index.bed` (325373)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 325373;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:18:47) elapsed: 7.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam`
1814382.0 704966

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam`
1814382.0 655479

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_gene_sort.bed` (329280,329284)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 329280;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 329284;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_gene_coverage.bed` (329856)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.204GB.  
  PID: 329856;	Command: bedtools;	Return code: 0;	Memory used: 0.031GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/hg38_annotations.bed.gz` (331694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331694;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/hg38_annotations.bed` (331697)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331697;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:18:58) elapsed: 10.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/raw/hg38_annotations.bed` (331714)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331714;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Enhancer_sort.bed` (331762,331763,331764,331765)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331762;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331763;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 331765;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB  
  PID: 331764;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Enhancer_plus_coverage.bed` (331783)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331783;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Enhancer_minus_coverage.bed` (331788)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331788;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_sort.bed` (331801,331802,331803,331804)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331801;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331802;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 331804;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 331803;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_plus_coverage.bed` (331807)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331807;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_minus_coverage.bed` (331816)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331816;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_Flanking_Region"` (331821)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331821;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed` (331822,331823,331824,331825)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331822;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331824;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 331823;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 331825;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_Flanking_Region_plus_coverage.bed` (331828)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331828;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_Flanking_Region_minus_coverage.bed` (331843)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331843;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/5_UTR"` (331853)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331853;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/5_UTR_sort.bed` (331854,331855,331856,331857)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331854;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331855;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 331857;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 331856;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_5_UTR_plus_coverage.bed` (331860)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331860;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_5_UTR_minus_coverage.bed` (331863)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331863;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/3_UTR"` (331868)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331868;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/3_UTR_sort.bed` (331869,331870,331871,331872)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331869;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331870;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 331872;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 331871;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_3_UTR_plus_coverage.bed` (331875)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331875;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_3_UTR_minus_coverage.bed` (331884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331884;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Exon_sort.bed` (331889,331890,331891,331892)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331889;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331890;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 331892;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 331891;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Exon_plus_coverage.bed` (331899)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331899;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Exon_minus_coverage.bed` (331910)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331910;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Intron_sort.bed` (331918,331919,331920,331921)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.204GB.  
  PID: 331918;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331920;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 331919;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 331921;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Intron_plus_coverage.bed` (331927)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331927;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Intron_minus_coverage.bed` (331945)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331945;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB


### Plot cFRiF/FRiF (06-14 21:19:36) elapsed: 38.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_10 -z 3099922541 -n 994670 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Intron_plus_coverage.bed` (331959)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.44GB.  
  PID: 331959;	Command: Rscript;	Return code: 0;	Memory used: 0.44GB

> `cFRiF`	QC_hg38/H9_PRO-seq_10_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_10_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_10 -z 3099922541 -n 994670 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_Intron_plus_coverage.bed` (332292)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 0.44GB.  
  PID: 332292;	Command: Rscript;	Return code: 0;	Memory used: 0.44GB

> `FRiF`	QC_hg38/H9_PRO-seq_10_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_10_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:20:41) elapsed: 65.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_exons_sort.bed` (332350,332351)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.44GB.  
  PID: 332351;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB  
  PID: 332350;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_introns_sort.bed` (332358,332359,332360)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.44GB.  
  PID: 332358;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 332360;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 332359;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exons_coverage.bed` (332366)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.44GB.  
  PID: 332366;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_introns_coverage.bed` (332371)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.44GB.  
  PID: 332371;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/1.814382)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exons_rpkm.bed` (332376,332377,332378)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.44GB.  
  PID: 332376;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 332378;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 332377;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/1.814382)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_introns_rpkm.bed` (332380,332381,332382)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.44GB.  
  PID: 332380;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 332382;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 332381;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exon_intron_ratios.bed` (332385,332386,332387)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.44GB.  
  PID: 332385;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 332387;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 332386;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exon_intron_ratios.bed --annotate` (332393)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.44GB.  
  PID: 332393;	Command: Rscript;	Return code: 0;	Memory used: 0.263GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_10_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_10_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/QC_hg38/H9_PRO-seq_10_exon_intron_ratios.bed` (332414)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.44GB.  
  PID: 332414;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (06-14 21:21:06) elapsed: 25.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam` (332421)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.44GB.  
  PID: 332421;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge --scale 1814382.0` (332424)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_10_plus_cuttrace_j1my62qp'
Processing with 1 cores...
Discarding 141 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000224v1', 'chrUn_GL000226v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000218v1']
Keeping 54 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270726v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270730v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270747v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrEBV']
Reduce step (merge files)...
Merging 54 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_plus_exact_body_0-mer.bw'
Merging 54 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:17:20. Running peak memory: 1.85GB.  
  PID: 332424;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.85GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam` (135565)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 1.85GB.  
  PID: 135565;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_minus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge --scale 1814382.0` (135567)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/aligned_hg38/H9_PRO-seq_10_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_10_minus_cuttrace_28dlyiwy'
Processing with 1 cores...
stdin is empty of data
Discarding 141 chunk(s) of reads: ['chrM', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1']
Keeping 54 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr9_KI270718v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_KI270743v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270750v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 54 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_minus_exact_body_0-mer.bw'
Merging 54 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_10/signal_hg38/H9_PRO-seq_10_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:17:01. Running peak memory: 2.477GB.  
  PID: 135567;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.477GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:36:56
*  Total elapsed time (all runs):  3:48:32
*         Peak memory (this run):  2.4766 GB
*        Pipeline completed time: 2020-06-14 21:55:30
