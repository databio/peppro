### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name HEK_ARF_PRO-seq --genome hg38 --input /project/shefflab/data/guertin/fastq/SRR8608070_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/SRR8608070_PE2.fastq.gz --protocol PRO --umi-len 8 --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj40-15c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/
*  Pipeline started at:   (02-18 11:08:58) elapsed: 1.0 _TIME_

### Version log:

*       Python version:  3.6.5
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.8.9
*        Pipeline hash:  d6fefacd05f25e266e8a162fa20e48a6f6c830d2
*      Pipeline branch:  * dev
*        Pipeline date:  2020-02-17 11:24:59 -0500

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `cutadapt`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `8`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data/guertin/fastq/SRR8608070_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data/guertin/fastq/SRR8608070_PE2.fastq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `16000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed/peppro/paper/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `HEK_ARF_PRO-seq`
*              `scale`:  `False`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `paired`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `8`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /project/shefflab/data/guertin/fastq/SRR8608070_PE1.fastq.gz
Local input file: /project/shefflab/data/guertin/fastq/SRR8608070_PE2.fastq.gz

> `File_mb`	1302.45	PEPPRO	_RES_

> `Read_type`	paired	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:08:59) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/SRR8608070_PE1.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R1.fastq.gz` (53173)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 53173;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/SRR8608070_PE2.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R2.fastq.gz` (53178)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 53178;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1.fastq`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R1.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1.fastq` (53184)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 0.003GB.  
  PID: 53184;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R2.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2.fastq` (428181)
<pre>
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 0.003GB.  
  PID: 428181;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	53049272	PEPPRO	_RES_

> `Fastq_reads`	53049272	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/HEK_ARF_PRO-seq_R2.fastq.gz']

### FASTQ processing:  (02-18 11:11:35) elapsed: 156.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R1_cutadapt.txt` (271802)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 0.281GB.  
  PID: 271802;	Command: cutadapt;	Return code: 0;	Memory used: 0.281GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq` (240938,241848)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.281GB.  
  PID: 240938;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 241848;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	16312692	PEPPRO	_RES_

> `Trim_loss_rate`	69.25	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq` (354766)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
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
Command completed. Elapsed time: 0:00:41. Running peak memory: 0.281GB.  
  PID: 354766;	Command: fastqc;	Return code: 0;	Memory used: 0.156GB

> `FastQC report r1`	fastqc/HEK_ARF_PRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_dedup.fastq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq` (129922)
<pre>
[INFO][0m 778059 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 1.029GB.  
  PID: 129922;	Command: seqkit;	Return code: 0;	Memory used: 1.029GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq` (219872,221238)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 1.029GB.  
  PID: 219872;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 221238;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	21262101.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	8863359.0	PEPPRO	_RES_

> `Duplicate_reads`	778059.0	PEPPRO	_RES_

> `Pct_reads_too_short`	33.4156	PEPPRO	_RES_

> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_R2_cutadapt.txt` (139898)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 1.029GB.  
  PID: 139898;	Command: cutadapt;	Return code: 0;	Memory used: 0.282GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq` (2616,2705)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 1.029GB.  
  PID: 2616;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 2705;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	17394524	PEPPRO	_RES_

> `Trim_loss_rate`	67.21	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq` (25575)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
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
Command completed. Elapsed time: 0:00:44. Running peak memory: 1.029GB.  
  PID: 25575;	Command: fastqc;	Return code: 0;	Memory used: 0.153GB

> `FastQC report r1`	fastqc/HEK_ARF_PRO-seq_R2_trimmed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.histogram`  

> `fastq_pair -t 47744344 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq` (25666)
<pre>
Left paired: 17619354		Right paired: 17619354
Left single: 41924		Right single: 1037248
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:13. Running peak memory: 2.442GB.  
  PID: 25666;	Command: fastq_pair;	Return code: 0;	Memory used: 2.442GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_noadap.fastq.paired.fq -o HEK_ARF_PRO-seq -d /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt` (25806)
<pre>
</pre>
Command completed. Elapsed time: 0:00:48. Running peak memory: 2.442GB.  
  PID: 25806;	Command: flash;	Return code: 0;	Memory used: 0.081GB


### Plot adapter insertion distribution (02-18 11:17:53) elapsed: 378.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt -u 8` (26016)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 2.442GB.  
  PID: 26016;	Command: Rscript;	Return code: 0;	Memory used: 0.205GB

> `Adapter insertion distribution`	cutadapt/HEK_ARF_PRO-seq_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/HEK_ARF_PRO-seq_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	17	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:18:12) elapsed: 19.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/cutadapt/HEK_ARF_PRO-seq.hist`

> `Degradation_ratio`	1.349	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq` (26063)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 2.442GB.  
  PID: 26063;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq.paired.fq`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq.paired.fq`  

> `fastq_pair -t 47744344 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq` (26088)
<pre>
Left paired: 16265754		Right paired: 16265754
Left single: 46939		Right single: 1128770
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 2.568GB.  
  PID: 26088;	Command: fastq_pair;	Return code: 0;	Memory used: 2.568GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/repaired.flag`  

> `mv /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq` (26214)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.568GB.  
  PID: 26214;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq` (26216)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.568GB.  
  PID: 26216;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/repaired.flag` (26217)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.568GB.  
  PID: 26217;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq.paired.fq`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq.paired.fq`  

> `fastq_pair -t 47744344 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq` (26218)
<pre>
Left paired: 15570095		Right paired: 15570095
Left single: 33539		Right single: 1824429
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 2.568GB.  
  PID: 26218;	Command: fastq_pair;	Return code: 0;	Memory used: 1.994GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/dups_repaired.flag`  

> `mv /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq` (26560)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 2.568GB.  
  PID: 26560;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq` (26561)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.568GB.  
  PID: 26561;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/dups_repaired.flag` (26562)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.568GB.  
  PID: 26562;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-18 11:20:38) elapsed: 146.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:20:38) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_bt2` (26563)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.568GB.  
  PID: 26563;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq` (26564)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HEK_ARF_PRO-seq -U /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
16265754 reads; of these:
  16265754 (100.00%) were unpaired; of these:
    14454678 (88.87%) aligned 0 times
    1811076 (11.13%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
11.13% overall alignment rate

> `Aligned_reads_human_rDNA`	3622152.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	20.82	PEPPRO	_RES_

### Map to human_rDNA (02-18 11:22:42) elapsed: 125.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_dups_bt2` (26802)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.568GB.  
  PID: 26802;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_dups_R2.fq` (26803)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HEK_ARF_PRO-seq -U /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/fastq/HEK_ARF_PRO-seq_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
1811076 reads skipped
0 reads lost

### Map to genome (02-18 11:24:29) elapsed: 107.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id HEK_ARF_PRO-seq -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/tmpaam2zkr0 -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam` (26968,26969,26970)
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
[bam_sort_core] merging from 6 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:39:37. Running peak memory: 3.654GB.  
  PID: 26968;	Command: bowtie2;	Return code: 0;	Memory used: 3.654GB  
  PID: 26969;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 26970;	Command: samtools;	Return code: 0;	Memory used: 0.88GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam` (32009)
<pre>
</pre>
Command completed. Elapsed time: 0:00:58. Running peak memory: 3.654GB.  
  PID: 32009;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	24911965	PEPPRO	_RES_

> `QC_filtered_reads`	15309473	PEPPRO	_RES_

> `Aligned_reads`	9602492.0	PEPPRO	_RES_

> `Alignment_rate`	55.2	PEPPRO	_RES_

> `Total_efficiency`	18.1	PEPPRO	_RES_

> `Read_depth`	2.05	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id HEK_ARF_PRO-seq -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/tmpaam2zkr0 -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp_dups.bam` (33314,33320,33324)
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
[bam_sort_core] merging from 6 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:40:50. Running peak memory: 3.654GB.  
  PID: 33320;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 33314;	Command: bowtie2;	Return code: 0;	Memory used: 3.653GB  
  PID: 33324;	Command: samtools;	Return code: 0;	Memory used: 0.887GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp_dups.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam` (37562)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 3.654GB.  
  PID: 37562;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (02-18 12:56:56) elapsed: 5547.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R2.fq` (37630)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.654GB.  
  PID: 37630;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/prealignments/HEK_ARF_PRO-seq_human_rDNA_unmap_R1.fq` (37656)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.654GB.  
  PID: 37656;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam` (37690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.654GB.  
  PID: 37690;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	252493	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam` (37711)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.654GB.  
  PID: 37711;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_noMT.bam` (37724,37725,37726,37727)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.654GB.  
  PID: 37726;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 37724;	Command: samtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 37725;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 37727;	Command: xargs;	Return code: 0;	Memory used: 0.054GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam` (37757)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 37757;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam` (37758)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.654GB.  
  PID: 37758;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Split BAM file (02-18 12:58:37) elapsed: 101.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam` (37792,37793)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:17. Running peak memory: 3.654GB.  
  PID: 37793;	Command: samtools;	Return code: 0;	Memory used: 2.039GB  
  PID: 37792;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE2.bam` (37908,37909)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:02. Running peak memory: 3.654GB.  
  PID: 37908;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 37909;	Command: samtools;	Return code: 0;	Memory used: 1.805GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp_dups.bam` (38269)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.654GB.  
  PID: 38269;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE1.bam` (38299,38300)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:17. Running peak memory: 3.654GB.  
  PID: 38300;	Command: samtools;	Return code: 0;	Memory used: 2.019GB  
  PID: 38299;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE2.bam` (38438,38439)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 3.654GB.  
  PID: 38438;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 38439;	Command: samtools;	Return code: 0;	Memory used: 1.687GB


### Calculate library complexity (02-18 13:03:56) elapsed: 319.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_out.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE1.bam` (38527)
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
Command completed. Elapsed time: 0:00:58. Running peak memory: 3.654GB.  
  PID: 38527;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE1.bam` (38587)
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
....................................................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:01. Running peak memory: 3.654GB.  
  PID: 38587;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam) > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_counts.txt` (38841)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.654GB.  
  PID: 38841;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_plot` (38861)
<pre>
Processing HEK_ARF_PRO-seq
INFO: Found real counts for HEK_ARF_PRO-seq - Total (M): 10.869678 Unique (M): 10.681803

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.654GB.  
  PID: 38861;	Command: Rscript;	Return code: 0;	Memory used: 0.205GB

> `Library complexity`	QC_hg38/HEK_ARF_PRO-seq_preseq_plot.pdf	Library complexity	QC_hg38/HEK_ARF_PRO-seq_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.9166	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 13:06:16) elapsed: 140.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam` (38881)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.654GB.  
  PID: 38881;	Command: samtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_bamQC.tsv` (38890)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/tmp_HEK_ARF_PRO-seq_PE1_yrts7lly'
Processing with 8 cores...
Discarding 96 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270757v1']
Keeping 99 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270508v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 3.654GB.  
  PID: 38890;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.82GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	5434839.0	PEPPRO	_RES_

> `PBC2`	5434839.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_unmap.bam` (38935)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.654GB.  
  PID: 38935;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_temp.bam`

> `Unmapped_reads`	3997391	PEPPRO	_RES_

### Split BAM by strand (02-18 13:06:51) elapsed: 35.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam` (38965)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.654GB.  
  PID: 38965;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam` (38998)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 3.654GB.  
  PID: 38998;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 13:07:58) elapsed: 68.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (39035)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39035;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_plus_TssEnrichment.txt` (39036)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.654GB.  
  PID: 39036;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.469GB


> `TSS_Plus_Score`	18.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_minus_TssEnrichment.txt` (39063)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.654GB.  
  PID: 39063;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.634GB


> `TSS_Minus_Score`	2.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_minus_TssEnrichment.txt` (39086)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.654GB.  
  PID: 39086;	Command: Rscript;	Return code: 0;	Memory used: 0.277GB

> `TSS enrichment`	QC_hg38/HEK_ARF_PRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/HEK_ARF_PRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt` (39106,39107,39108,39109)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39106;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 39108;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 39107;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 39109;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt` (39111)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39111;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (02-18 13:08:14) elapsed: 16.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_tss.bed` (39113,39114)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.654GB.  
  PID: 39113;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 39114;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (39117,39118)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39117;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 39118;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_TSS_density.bed` (39121,39122,39123,39124)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.654GB.  
  PID: 39124;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 39121;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 39123;	Command: sort;	Return code: 0;	Memory used: 0.008GB  
  PID: 39122;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_gene_body_density.bed` (39144,39145,39146)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.654GB.  
  PID: 39144;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB  
  PID: 39146;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 39145;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed` (39166,39167,39168)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39166;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 39168;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 39167;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	17.17	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed` (39174)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.654GB.  
  PID: 39174;	Command: Rscript;	Return code: 0;	Memory used: 0.299GB

> `Pause index`	QC_hg38/HEK_ARF_PRO-seq_pause_index.pdf	Pause index	QC_hg38/HEK_ARF_PRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_pause_index.bed` (39192)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39192;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 13:08:56) elapsed: 42.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam`
9602492.0 4000690

> `Plus_FRiP`	0.42	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam`
9602492.0 3788384

> `Minus_FRiP`	0.39	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_gene_sort.bed` (39234,39235)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.654GB.  
  PID: 39234;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 39235;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_gene_coverage.bed` (39237)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.654GB.  
  PID: 39237;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed.gz` (39262)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39262;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed` (39263)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39263;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 13:09:37) elapsed: 41.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/raw/hg38_annotations.bed` (39273)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.654GB.  
  PID: 39273;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR"` (39275)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39275;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR_sort.bed` (39276,39278,39279,39280)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.654GB.  
  PID: 39276;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 39278;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 39280;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 39279;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_plus_coverage.bed` (39282)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.654GB.  
  PID: 39282;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_minus_coverage.bed` (39291)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.654GB.  
  PID: 39291;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR"` (39298)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39298;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR_sort.bed` (39299,39300,39301,39302)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.654GB.  
  PID: 39299;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 39300;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 39302;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 39301;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_plus_coverage.bed` (39305)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.654GB.  
  PID: 39305;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_minus_coverage.bed` (39312)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.654GB.  
  PID: 39312;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer_sort.bed` (39537,39538,39539,39540)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.654GB.  
  PID: 39537;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 39538;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 39540;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB  
  PID: 39539;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_plus_coverage.bed` (39543)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.654GB.  
  PID: 39543;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_minus_coverage.bed` (39554)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.654GB.  
  PID: 39554;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon_sort.bed` (39564,39565,39566,39567)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.654GB.  
  PID: 39564;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 39565;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 39567;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 39566;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_plus_coverage.bed` (39573)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.654GB.  
  PID: 39573;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_minus_coverage.bed` (39587)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.654GB.  
  PID: 39587;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron_sort.bed` (39598,39599,39600,39601)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.654GB.  
  PID: 39598;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 39600;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 39599;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 39601;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_plus_coverage.bed` (39605)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.654GB.  
  PID: 39605;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_minus_coverage.bed` (39614)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.654GB.  
  PID: 39614;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_sort.bed` (39623,39624,39625,39626)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39623;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 39624;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 39626;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 39625;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_plus_coverage.bed` (39628)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.654GB.  
  PID: 39628;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_minus_coverage.bed` (39639)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.654GB.  
  PID: 39639;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region"` (39646)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39646;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (39647,39649,39650,39651)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.654GB.  
  PID: 39647;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 39650;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 39649;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 39651;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (39653)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.654GB.  
  PID: 39653;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (39661)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.654GB.  
  PID: 39661;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


### Plot FRiF/PRiF (02-18 13:11:30) elapsed: 113.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s HEK_ARF_PRO-seq -z 3099922541 -n 5523870 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (39685)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.654GB.  
  PID: 39685;	Command: Rscript;	Return code: 0;	Memory used: 0.469GB

> `FRiF`	QC_hg38/HEK_ARF_PRO-seq_frif.pdf	FRiF	QC_hg38/HEK_ARF_PRO-seq_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s HEK_ARF_PRO-seq -z 3099922541 -n 5523870 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (39732)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.654GB.  
  PID: 39732;	Command: Rscript;	Return code: 0;	Memory used: 0.469GB

> `PRiF`	QC_hg38/HEK_ARF_PRO-seq_prif.pdf	PRiF	QC_hg38/HEK_ARF_PRO-seq_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 13:12:33) elapsed: 63.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_exons_sort.bed` (39772,39773)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.654GB.  
  PID: 39773;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB  
  PID: 39772;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_introns_sort.bed` (39780,39781,39782)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.654GB.  
  PID: 39780;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 39782;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB  
  PID: 39781;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_coverage.bed` (39789)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 3.654GB.  
  PID: 39789;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_coverage.bed` (39802)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.654GB.  
  PID: 39802;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/9.602492)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_rpkm.bed` (39817,39819,39820)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.654GB.  
  PID: 39817;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 39820;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 39819;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/9.602492)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_rpkm.bed` (39822,39823,39824)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.654GB.  
  PID: 39822;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 39824;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 39823;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed` (39827,39828,39829)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39827;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 39829;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 39828;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed --annotate` (39835)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.654GB.  
  PID: 39835;	Command: Rscript;	Return code: 0;	Memory used: 0.238GB

> `mRNA contamination`	QC_hg38/HEK_ARF_PRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/HEK_ARF_PRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/QC_hg38/HEK_ARF_PRO-seq_exon_intron_ratios.bed` (39857)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.654GB.  
  PID: 39857;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-18 13:13:20) elapsed: 47.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam` (39865)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.654GB.  
  PID: 39865;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (39870)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_HEK_ARF_PRO-seq_plus_cuttrace_mw_r_gie'
Processing with 2 cores...
Discarding 113 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270752v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 82 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270589v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_exact_body_0-mer.bw'
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:08. Running peak memory: 3.654GB.  
  PID: 39870;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.053GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam` (41674)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.654GB.  
  PID: 41674;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_minus_exact_body_0-mer.bw -p 5 --variable-step --tail-edge` (41679)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/aligned_hg38/HEK_ARF_PRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_HEK_ARF_PRO-seq_minus_cuttrace_olyfi1zo'
Processing with 5 cores...
stdin is empty of data
Discarding 107 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_KI270722v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 88 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270508v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 88 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/HEK_ARF_PRO-seq/signal_hg38/HEK_ARF_PRO-seq_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.654GB.  
  PID: 41679;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.235GB

Starting cleanup: 68 files; 4 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:14:51
*  Total elapsed time (all runs):  4:36:34
*         Peak memory (this run):  3.6543 GB
*        Pipeline completed time: 2020-02-18 13:23:49
