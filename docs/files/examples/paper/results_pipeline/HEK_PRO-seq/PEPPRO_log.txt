### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name HEK_PRO-seq --genome hg38 --input /project/shefflab/data//guertin/fastq/SRR8608074_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/SRR8608074_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-32c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/
*  Pipeline started at:   (06-15 07:17:12) elapsed: 0.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/SRR8608074_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/SRR8608074_PE2.fastq.gz']`
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
*        `sample_name`:  `HEK_PRO-seq`
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

Local input file: /project/shefflab/data//guertin/fastq/SRR8608074_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/SRR8608074_PE2.fastq.gz

> `File_mb`	1337.94	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:13) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/SRR8608074_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R1.fastq.gz` (380687)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 380687;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/SRR8608074_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R2.fastq.gz` (380691)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 380691;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1.fastq` (380697)
<pre>
</pre>
Command completed. Elapsed time: 0:01:45. Running peak memory: 0.002GB.  
  PID: 380697;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2.fastq` (380803)
<pre>
</pre>
Command completed. Elapsed time: 0:00:59. Running peak memory: 0.002GB.  
  PID: 380803;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	51533198	PEPPRO	_RES_

> `Fastq_reads`	51533198	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R2.fastq.gz']

### FASTQ processing:  (06-15 07:20:20) elapsed: 187.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R1_cutadapt.txt` (381203)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 1.773GB.  
  PID: 381203;	Command: cutadapt;	Return code: 0;	Memory used: 1.773GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq` (381274,381275)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 1.773GB.  
  PID: 381274;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 381275;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	20103582	PEPPRO	_RES_

> `Trim_loss_rate`	60.99	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq` (381328)
<pre>
Started analysis of HEK_PRO-seq_R1_processed.fastq
Approx 5% complete for HEK_PRO-seq_R1_processed.fastq
Approx 10% complete for HEK_PRO-seq_R1_processed.fastq
Approx 15% complete for HEK_PRO-seq_R1_processed.fastq
Approx 20% complete for HEK_PRO-seq_R1_processed.fastq
Approx 25% complete for HEK_PRO-seq_R1_processed.fastq
Approx 30% complete for HEK_PRO-seq_R1_processed.fastq
Approx 35% complete for HEK_PRO-seq_R1_processed.fastq
Approx 40% complete for HEK_PRO-seq_R1_processed.fastq
Approx 45% complete for HEK_PRO-seq_R1_processed.fastq
Approx 50% complete for HEK_PRO-seq_R1_processed.fastq
Approx 55% complete for HEK_PRO-seq_R1_processed.fastq
Approx 60% complete for HEK_PRO-seq_R1_processed.fastq
Approx 65% complete for HEK_PRO-seq_R1_processed.fastq
Approx 70% complete for HEK_PRO-seq_R1_processed.fastq
Approx 75% complete for HEK_PRO-seq_R1_processed.fastq
Approx 80% complete for HEK_PRO-seq_R1_processed.fastq
Approx 85% complete for HEK_PRO-seq_R1_processed.fastq
Approx 90% complete for HEK_PRO-seq_R1_processed.fastq
Approx 95% complete for HEK_PRO-seq_R1_processed.fastq
Analysis complete for HEK_PRO-seq_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 1.773GB.  
  PID: 381328;	Command: fastqc;	Return code: 0;	Memory used: 0.156GB

> `FastQC report r1`	fastqc/HEK_PRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq` (381432)
<pre>
[INFO][0m 916082 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 1.986GB.  
  PID: 381432;	Command: seqkit;	Return code: 0;	Memory used: 1.986GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq` (381515,381516)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 1.986GB.  
  PID: 381515;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 381516;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	16626152.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	5380555.0	PEPPRO	_RES_

> `Duplicate_reads`	916082.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	20.8819	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/processed_R1.flag` (381619)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.986GB.  
  PID: 381619;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R2_cutadapt.txt` (381621)
<pre>
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 1.986GB.  
  PID: 381621;	Command: cutadapt;	Return code: 0;	Memory used: 1.927GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq` (381919,381920)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 1.986GB.  
  PID: 381919;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 381920;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	40207164	PEPPRO	_RES_

> `Trim_loss_rate`	21.98	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq` (381983)
<pre>
Started analysis of HEK_PRO-seq_R1_processed.fastq
Approx 5% complete for HEK_PRO-seq_R1_processed.fastq
Approx 10% complete for HEK_PRO-seq_R1_processed.fastq
Approx 15% complete for HEK_PRO-seq_R1_processed.fastq
Approx 20% complete for HEK_PRO-seq_R1_processed.fastq
Approx 25% complete for HEK_PRO-seq_R1_processed.fastq
Approx 30% complete for HEK_PRO-seq_R1_processed.fastq
Approx 35% complete for HEK_PRO-seq_R1_processed.fastq
Approx 40% complete for HEK_PRO-seq_R1_processed.fastq
Approx 45% complete for HEK_PRO-seq_R1_processed.fastq
Approx 50% complete for HEK_PRO-seq_R1_processed.fastq
Approx 55% complete for HEK_PRO-seq_R1_processed.fastq
Approx 60% complete for HEK_PRO-seq_R1_processed.fastq
Approx 65% complete for HEK_PRO-seq_R1_processed.fastq
Approx 70% complete for HEK_PRO-seq_R1_processed.fastq
Approx 75% complete for HEK_PRO-seq_R1_processed.fastq
Approx 80% complete for HEK_PRO-seq_R1_processed.fastq
Approx 85% complete for HEK_PRO-seq_R1_processed.fastq
Approx 90% complete for HEK_PRO-seq_R1_processed.fastq
Approx 95% complete for HEK_PRO-seq_R1_processed.fastq
Analysis complete for HEK_PRO-seq_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 1.986GB.  
  PID: 381983;	Command: fastqc;	Return code: 0;	Memory used: 0.155GB

> `FastQC report r1`	fastqc/HEK_PRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq` (382054)
<pre>
Started analysis of HEK_PRO-seq_R2_trimmed.fastq
Approx 5% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 10% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 15% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 20% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 25% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 30% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 35% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 40% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 45% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 50% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 55% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 60% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 65% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 70% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 75% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 80% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 85% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 90% complete for HEK_PRO-seq_R2_trimmed.fastq
Approx 95% complete for HEK_PRO-seq_R2_trimmed.fastq
Analysis complete for HEK_PRO-seq_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 1.986GB.  
  PID: 382054;	Command: fastqc;	Return code: 0;	Memory used: 0.153GB

> `FastQC report r2`	fastqc/HEK_PRO-seq_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.histogram`  

> `fastq_pair -t 46379878 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq` (382246)
<pre>
Left paired: 20354925		Right paired: 20354925
Left single: 31120		Right single: 456238
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 2.441GB.  
  PID: 382246;	Command: fastq_pair;	Return code: 0;	Memory used: 2.441GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq.paired.fq -o HEK_PRO-seq -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt` (382390)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 2.441GB.  
  PID: 382390;	Command: flash;	Return code: 0;	Memory used: 0.133GB


### Plot adapter insertion distribution (06-15 07:29:38) elapsed: 558.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt -u 8` (382583)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 2.441GB.  
  PID: 382583;	Command: Rscript;	Return code: 0;	Memory used: 0.287GB

> `Adapter insertion distribution`	cutadapt/HEK_PRO-seq_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/HEK_PRO-seq_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	37	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:29:44) elapsed: 6.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`

> `Degradation_ratio`	0.5783	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq` (382613)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 2.441GB.  
  PID: 382613;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/processed_R2.flag` (382911)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.441GB.  
  PID: 382911;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/repaired.flag`  

> `fastq_pair -t 46379878 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq` (382912)
<pre>
Left paired: 20071296		Right paired: 20071296
Left single: 32287		Right single: 474552
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:46. Running peak memory: 2.441GB.  
  PID: 382912;	Command: fastq_pair;	Return code: 0;	Memory used: 2.009GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq` (383044)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.441GB.  
  PID: 383044;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq` (383045)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.441GB.  
  PID: 383045;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/repaired.flag` (383047)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.441GB.  
  PID: 383047;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/dups_repaired.flag`  

> `fastq_pair -t 46379878 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq` (383048)
<pre>
Left paired: 19171561		Right paired: 19171561
Left single: 25508		Right single: 1374287
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:18. Running peak memory: 2.593GB.  
  PID: 383048;	Command: fastq_pair;	Return code: 0;	Memory used: 2.593GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq` (383175)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 2.593GB.  
  PID: 383175;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq` (383182)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 2.593GB.  
  PID: 383182;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/dups_repaired.flag` (383183)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.593GB.  
  PID: 383183;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:33:13) elapsed: 209.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:33:13) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_bt2` (383184)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.593GB.  
  PID: 383184;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq` (383186)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HEK_PRO-seq -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
20071296 reads; of these:
  20071296 (100.00%) were unpaired; of these:
    17956185 (89.46%) aligned 0 times
    2115111 (10.54%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.54% overall alignment rate

> `Aligned_reads_human_rDNA`	4230222.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.52	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:36:01) elapsed: 168.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_dups_bt2` (383576)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.593GB.  
  PID: 383576;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_dups_R2.fq` (383578)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HEK_PRO-seq -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2115112 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:38:40) elapsed: 159.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id HEK_PRO-seq -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/tmp6ifjo32x -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam` (383744,383745,383746)
<pre>
1726816 reads skipped
0 reads lost
17956184 reads; of these:
  17956184 (100.00%) were paired; of these:
    5503120 (30.65%) aligned concordantly 0 times
    10297215 (57.35%) aligned concordantly exactly 1 time
    2155849 (12.01%) aligned concordantly >1 times
    ----
    5503120 pairs aligned concordantly 0 times; of these:
      1505429 (27.36%) aligned discordantly 1 time
    ----
    3997691 pairs aligned 0 times concordantly or discordantly; of these:
      7995382 mates make up the pairs; of these:
        2807990 (35.12%) aligned 0 times
        2063628 (25.81%) aligned exactly 1 time
        3123764 (39.07%) aligned >1 times
92.18% overall alignment rate
[bam_sort_core] merging from 9 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:33:26. Running peak memory: 3.772GB.  
  PID: 383745;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 383744;	Command: bowtie2;	Return code: 0;	Memory used: 3.772GB  
  PID: 383746;	Command: samtools;	Return code: 0;	Memory used: 0.928GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam` (387133)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 3.772GB.  
  PID: 387133;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	33104378	PEPPRO	_RES_

> `QC_filtered_reads`	18769818	PEPPRO	_RES_

> `Aligned_reads`	14334559.5	PEPPRO	_RES_

> `Alignment_rate`	35.65	PEPPRO	_RES_

> `Total_efficiency`	27.82	PEPPRO	_RES_

> `Read_depth`	2.41	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id HEK_PRO-seq -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/tmp6ifjo32x -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp_dups.bam` (388990,388996,388997)
<pre>
17444745 reads; of these:
  17444745 (100.00%) were paired; of these:
    5326191 (30.53%) aligned concordantly 0 times
    10028868 (57.49%) aligned concordantly exactly 1 time
    2089686 (11.98%) aligned concordantly >1 times
    ----
    5326191 pairs aligned concordantly 0 times; of these:
      1467659 (27.56%) aligned discordantly 1 time
    ----
    3858532 pairs aligned 0 times concordantly or discordantly; of these:
      7717064 mates make up the pairs; of these:
        2726098 (35.33%) aligned 0 times
        1998944 (25.90%) aligned exactly 1 time
        2992022 (38.77%) aligned >1 times
92.19% overall alignment rate
[bam_sort_core] merging from 9 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:30:01. Running peak memory: 3.772GB.  
  PID: 388990;	Command: bowtie2;	Return code: 0;	Memory used: 3.756GB  
  PID: 388996;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 388997;	Command: samtools;	Return code: 0;	Memory used: 0.906GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam` (393207)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 3.772GB.  
  PID: 393207;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 08:57:19) elapsed: 4719.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R1.fq` (393916)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.772GB.  
  PID: 393916;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq` (393975)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.772GB.  
  PID: 393975;	Command: pigz;	Return code: 0;	Memory used: 0.013GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam` (394034)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.772GB.  
  PID: 394034;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	325556	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam` (394095)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.772GB.  
  PID: 394095;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/chr_sizes.bed` (394141,394142,394143,394144)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 394142;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 394144;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 394141;	Command: samtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 394143;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_noMT.bam` (394146)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 3.772GB.  
  PID: 394146;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam` (394211)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 394211;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam` (394213)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.772GB.  
  PID: 394213;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Split BAM file (06-15 08:59:34) elapsed: 135.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam` (394562,394563)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:44. Running peak memory: 3.772GB.  
  PID: 394562;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 394563;	Command: samtools;	Return code: 0;	Memory used: 3.641GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE2.bam` (395221,395222)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:39. Running peak memory: 3.772GB.  
  PID: 395221;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 395222;	Command: samtools;	Return code: 0;	Memory used: 3.076GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp_dups.bam` (395571)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.772GB.  
  PID: 395571;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE1.bam` (396181,396182)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:53. Running peak memory: 3.772GB.  
  PID: 396181;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 396182;	Command: samtools;	Return code: 0;	Memory used: 3.58GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE2.bam` (396694,396695)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:38. Running peak memory: 3.772GB.  
  PID: 396694;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 396695;	Command: samtools;	Return code: 0;	Memory used: 3.022GB


### Calculate library complexity (06-15 09:07:34) elapsed: 479.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE1.bam` (397460)
<pre>
BAM_INPUT
TOTAL READS     = 14923343
COUNTS_SUM      = 14923343
DISTINCT READS  = 1.32013e+07
DISTINCT COUNTS = 220
MAX COUNT       = 14762
COUNTS OF 1     = 1.216e+07
OBSERVED COUNTS (14763)
1	12160007
2	801020
3	141654
4	44209
5	19541
6	10184
7	6080
8	3969
9	2710
10	1970
11	1530
12	1137
13	864
14	760
15	571
16	491
17	498
18	370
19	342
20	262
21	230
22	212
23	195
24	181
25	170
26	146
27	128
28	114
29	96
30	109
31	85
32	81
33	65
34	71
35	68
36	72
37	52
38	39
39	48
40	40
41	41
42	27
43	34
44	29
45	27
46	38
47	28
48	36
49	36
50	21
51	18
52	24
53	12
54	23
55	15
56	18
57	21
58	15
59	18
60	8
61	14
62	16
63	11
64	8
65	8
66	10
67	10
68	11
69	12
70	12
71	5
72	6
73	10
74	5
75	8
76	10
77	5
78	8
79	8
80	5
81	7
82	4
83	5
84	6
85	6
86	3
87	5
88	5
89	4
90	2
91	4
92	6
93	4
94	6
95	5
96	9
97	3
98	6
99	2
100	2
101	4
103	3
104	2
105	4
106	1
107	2
108	3
109	2
110	5
111	4
112	7
113	2
114	3
115	2
116	1
117	1
118	1
120	2
121	2
122	4
123	1
124	1
125	1
126	5
127	2
128	1
129	1
130	3
132	4
133	4
134	2
135	1
137	2
138	1
139	3
140	2
141	2
142	2
145	1
147	2
148	5
149	1
152	2
153	1
158	1
159	1
160	1
161	2
163	1
164	1
165	1
166	1
168	2
169	1
170	1
172	2
176	1
177	1
179	1
182	1
184	2
185	2
188	1
189	1
191	3
194	1
202	1
206	2
208	1
209	2
211	1
215	2
217	1
222	1
224	1
225	3
226	2
229	1
240	2
241	1
244	1
248	1
275	1
276	1
277	1
283	1
287	1
288	1
289	1
292	1
328	1
351	1
358	1
361	1
363	1
404	1
431	2
436	1
469	1
544	1
581	1
610	1
612	1
696	1
726	1
758	1
778	1
880	1
973	1
1110	1
1493	1
1596	1
1700	1
1764	1
1771	1
2430	1
3752	1
5609	1
8568	1
14762	1

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
</pre>
Command completed. Elapsed time: 0:01:19. Running peak memory: 3.772GB.  
  PID: 397460;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE1.bam` (397570)
<pre>
BAM_INPUT
TOTAL READS     = 14923343
DISTINCT READS  = 1.32013e+07
DISTINCT COUNTS = 220
MAX COUNT       = 14762
COUNTS OF 1     = 1.216e+07
MAX TERMS       = 100
OBSERVED COUNTS (14763)
1	12160007
2	801020
3	141654
4	44209
5	19541
6	10184
7	6080
8	3969
9	2710
10	1970
11	1530
12	1137
13	864
14	760
15	571
16	491
17	498
18	370
19	342
20	262
21	230
22	212
23	195
24	181
25	170
26	146
27	128
28	114
29	96
30	109
31	85
32	81
33	65
34	71
35	68
36	72
37	52
38	39
39	48
40	40
41	41
42	27
43	34
44	29
45	27
46	38
47	28
48	36
49	36
50	21
51	18
52	24
53	12
54	23
55	15
56	18
57	21
58	15
59	18
60	8
61	14
62	16
63	11
64	8
65	8
66	10
67	10
68	11
69	12
70	12
71	5
72	6
73	10
74	5
75	8
76	10
77	5
78	8
79	8
80	5
81	7
82	4
83	5
84	6
85	6
86	3
87	5
88	5
89	4
90	2
91	4
92	6
93	4
94	6
95	5
96	9
97	3
98	6
99	2
100	2
101	4
103	3
104	2
105	4
106	1
107	2
108	3
109	2
110	5
111	4
112	7
113	2
114	3
115	2
116	1
117	1
118	1
120	2
121	2
122	4
123	1
124	1
125	1
126	5
127	2
128	1
129	1
130	3
132	4
133	4
134	2
135	1
137	2
138	1
139	3
140	2
141	2
142	2
145	1
147	2
148	5
149	1
152	2
153	1
158	1
159	1
160	1
161	2
163	1
164	1
165	1
166	1
168	2
169	1
170	1
172	2
176	1
177	1
179	1
182	1
184	2
185	2
188	1
189	1
191	3
194	1
202	1
206	2
208	1
209	2
211	1
215	2
217	1
222	1
224	1
225	3
226	2
229	1
240	2
241	1
244	1
248	1
275	1
276	1
277	1
283	1
287	1
288	1
289	1
292	1
328	1
351	1
358	1
361	1
363	1
404	1
431	2
436	1
469	1
544	1
581	1
610	1
612	1
696	1
726	1
758	1
778	1
880	1
973	1
1110	1
1493	1
1596	1
1700	1
1764	1
1771	1
2430	1
3752	1
5609	1
8568	1
14762	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
....................................................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:24. Running peak memory: 3.772GB.  
  PID: 397570;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_counts.txt` (397921)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.772GB.  
  PID: 397921;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_plot` (397941)
<pre>
Processing HEK_PRO-seq
INFO: Found real counts for HEK_PRO-seq - Total (M): 15.183469 Unique (M): 14.923343

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.772GB.  
  PID: 397941;	Command: Rscript;	Return code: 0;	Memory used: 0.268GB

> `Library complexity`	QC_hg38/HEK_PRO-seq_preseq_plot.pdf	Library complexity	QC_hg38/HEK_PRO-seq_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.9099	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 09:10:42) elapsed: 188.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam` (397960)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.772GB.  
  PID: 397960;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_bamQC.tsv` (397979)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/tmp_HEK_PRO-seq_PE1_q52j1sna'
Processing with 12 cores...
Discarding 105 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270748v1', 'chrUn_KI270752v1', 'chrUn_KI270756v1']
Keeping 90 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.772GB.  
  PID: 397979;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.071GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	7591734.5	PEPPRO	_RES_

> `PBC2`	7591734.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_unmap.bam` (398020)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.772GB.  
  PID: 398020;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam`

> `Unmapped_reads`	2807990	PEPPRO	_RES_

### Split BAM by strand (06-15 09:11:17) elapsed: 35.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam` (398059)
<pre>
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 3.772GB.  
  PID: 398059;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam` (398099)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 3.772GB.  
  PID: 398099;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 09:12:46) elapsed: 89.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (398152)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398152;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_plus_TssEnrichment.txt` (398154)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.772GB.  
  PID: 398154;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.987GB


> `TSS_coding_score`	16.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_minus_TssEnrichment.txt` (398185)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.772GB.  
  PID: 398185;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.981GB


> `TSS_non-coding_score`	2.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_minus_TssEnrichment.txt` (398216)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.772GB.  
  PID: 398216;	Command: Rscript;	Return code: 0;	Memory used: 0.32GB

> `TSS enrichment`	QC_hg38/HEK_PRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/HEK_PRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt` (398238,398239,398240,398241)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398238;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 398240;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 398239;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 398241;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt` (398243)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398243;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 09:13:03) elapsed: 17.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_tss.bed` (398245,398246)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.772GB.  
  PID: 398245;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 398246;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (398250,398251)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398250;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 398251;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_TSS_density.bed` (398254,398255,398256,398257)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.772GB.  
  PID: 398254;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 398256;	Command: sort;	Return code: 0;	Memory used: 0.012GB  
  PID: 398255;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 398257;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_gene_body_density.bed` (398274,398275,398276)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.772GB.  
  PID: 398274;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 398276;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 398275;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/tmp3fy58s0e` (398295,398296,398297)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398295;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 398297;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 398296;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/tmp3fy58s0e | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0136883) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/tmp3fy58s0e > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed` (398304)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398304;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.23	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed` (398309)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.772GB.  
  PID: 398309;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/HEK_PRO-seq_pause_index.pdf	Pause index	QC_hg38/HEK_PRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed` (398330)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398330;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 09:13:50) elapsed: 47.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam`
14334559.5 5618447

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam`
14334559.5 5337871

> `Minus_FRiP`	0.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_gene_sort.bed` (398370,398371)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.772GB.  
  PID: 398370;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 398371;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_gene_coverage.bed` (398373)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.772GB.  
  PID: 398373;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed.gz` (398391)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398391;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed` (398392)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398392;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 09:14:38) elapsed: 49.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed` (398401)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.772GB.  
  PID: 398401;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer_sort.bed` (398403,398404,398405,398406)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.772GB.  
  PID: 398403;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 398404;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 398406;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB  
  PID: 398405;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_plus_coverage.bed` (398409)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.772GB.  
  PID: 398409;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_minus_coverage.bed` (398417)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.772GB.  
  PID: 398417;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_sort.bed` (398428,398429,398430,398431)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398428;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 398429;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 398431;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 398430;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_plus_coverage.bed` (398433)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.772GB.  
  PID: 398433;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_minus_coverage.bed` (398634)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.772GB.  
  PID: 398634;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region"` (398646)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398646;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (398647,398648,398649,398650)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.772GB.  
  PID: 398647;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 398649;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 398648;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 398650;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (398653)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.772GB.  
  PID: 398653;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (398662)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.772GB.  
  PID: 398662;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR"` (398671)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398671;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR_sort.bed` (398672,398673,398674,398675)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.772GB.  
  PID: 398672;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 398673;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 398675;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 398674;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_plus_coverage.bed` (398678)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.772GB.  
  PID: 398678;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_minus_coverage.bed` (398687)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.772GB.  
  PID: 398687;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR"` (398696)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398696;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR_sort.bed` (398697,398698,398699,398700)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.772GB.  
  PID: 398697;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 398698;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 398700;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 398699;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_plus_coverage.bed` (398702)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.772GB.  
  PID: 398702;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_minus_coverage.bed` (398711)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.772GB.  
  PID: 398711;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Exon_sort.bed` (398721,398722,398723,398724)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.772GB.  
  PID: 398721;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 398722;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 398724;	Command: bedtools;	Return code: 0;	Memory used: 0.17GB  
  PID: 398723;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_plus_coverage.bed` (398728)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.772GB.  
  PID: 398728;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_minus_coverage.bed` (398738)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.772GB.  
  PID: 398738;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Intron_sort.bed` (398748,398749,398750,398751)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.772GB.  
  PID: 398748;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 398750;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 398749;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 398751;	Command: bedtools;	Return code: 0;	Memory used: 0.082GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_plus_coverage.bed` (398754)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.772GB.  
  PID: 398754;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_minus_coverage.bed` (398766)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.772GB.  
  PID: 398766;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB


### Plot cFRiF/FRiF (06-15 09:16:55) elapsed: 137.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s HEK_PRO-seq -z 3099922541 -n 7704166 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_plus_coverage.bed` (398791)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.772GB.  
  PID: 398791;	Command: Rscript;	Return code: 0;	Memory used: 0.502GB

> `cFRiF`	QC_hg38/HEK_PRO-seq_cFRiF.pdf	cFRiF	QC_hg38/HEK_PRO-seq_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s HEK_PRO-seq -z 3099922541 -n 7704166 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_plus_coverage.bed` (398831)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.772GB.  
  PID: 398831;	Command: Rscript;	Return code: 0;	Memory used: 0.472GB

> `FRiF`	QC_hg38/HEK_PRO-seq_FRiF.pdf	FRiF	QC_hg38/HEK_PRO-seq_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 09:17:57) elapsed: 62.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_exons_sort.bed` (398864,398865)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.772GB.  
  PID: 398864;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 398865;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_introns_sort.bed` (398871,398872,398873)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.772GB.  
  PID: 398871;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 398873;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 398872;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_coverage.bed` (398880)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.772GB.  
  PID: 398880;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_coverage.bed` (398896)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.772GB.  
  PID: 398896;	Command: bedtools;	Return code: 0;	Memory used: 0.026GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.3345595)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_rpkm.bed` (398914,398915,398916)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.772GB.  
  PID: 398914;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 398916;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 398915;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.3345595)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_rpkm.bed` (398919,398920,398921)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.772GB.  
  PID: 398919;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 398921;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 398920;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed` (398923,398924,398925)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398923;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 398925;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 398924;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.25	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed --annotate` (398931)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.772GB.  
  PID: 398931;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/HEK_PRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/HEK_PRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed` (398952)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.772GB.  
  PID: 398952;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 09:18:49) elapsed: 52.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam` (398961)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.772GB.  
  PID: 398961;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 14334559.5` (398966)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_HEK_PRO-seq_plus_cuttrace_9fiv4kcf'
Processing with 4 cores...
Discarding 121 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 74 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 74 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_exact_body_0-mer.bw'
Merging 74 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:15. Running peak memory: 3.772GB.  
  PID: 398966;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.405GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam` (400419)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.772GB.  
  PID: 400419;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 14334559.5` (400428)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_HEK_PRO-seq_minus_cuttrace_5prlm0ru'
Processing with 4 cores...
Discarding 117 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr5_GL000208v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr17_KI270730v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 78 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 78 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_exact_body_0-mer.bw'
Merging 78 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:06:16. Running peak memory: 3.772GB.  
  PID: 400428;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.405GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:14:20
*  Total elapsed time (all runs):  3:55:34
*         Peak memory (this run):  3.7719 GB
*        Pipeline completed time: 2020-06-15 09:31:31
