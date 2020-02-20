### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name H9_PRO-seq_1 --genome hg38 --input /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz --protocol PRO --umi-len 8 --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj40-15c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/
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
*              `input`:  `['/project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_1`
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

Local input file: /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz
Local input file: /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz

> `File_mb`	2469.73	PEPPRO	_RES_

> `Read_type`	paired	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:08:59) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE1.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz` (53174)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 53174;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep1_PE2.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz` (53177)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 53177;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq` (53182)
<pre>
</pre>
Command completed. Elapsed time: 0:02:48. Running peak memory: 0.003GB.  
  PID: 53182;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2.fastq` (329689)
<pre>
</pre>
Command completed. Elapsed time: 0:01:43. Running peak memory: 0.003GB.  
  PID: 329689;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	97729860	PEPPRO	_RES_

> `Fastq_reads`	97729860	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/H9_PRO-seq_1_R2.fastq.gz']

### FASTQ processing:  (02-18 11:14:29) elapsed: 330.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt` (305296)
<pre>
</pre>
Command completed. Elapsed time: 0:01:50. Running peak memory: 0.281GB.  
  PID: 305296;	Command: cutadapt;	Return code: 0;	Memory used: 0.281GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (25735,25736)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 0.281GB.  
  PID: 25735;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 25736;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	30221343	PEPPRO	_RES_

> `Trim_loss_rate`	69.08	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (25784)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_1_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_1_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_1_R1_processed.fastq
Analysis complete for H9_PRO-seq_1_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:17. Running peak memory: 0.281GB.  
  PID: 25784;	Command: fastqc;	Return code: 0;	Memory used: 0.175GB

> `FastQC report r1`	fastqc/H9_PRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_dedup.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq` (26074)
<pre>
[INFO][0m 3516802 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 1.994GB.  
  PID: 26074;	Command: seqkit;	Return code: 0;	Memory used: 1.994GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq` (26203,26204)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 1.994GB.  
  PID: 26203;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 26204;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	34822522.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	17998987.0	PEPPRO	_RES_

> `Duplicate_reads`	3516802.0	PEPPRO	_RES_

> `Pct_reads_too_short`	36.8342	PEPPRO	_RES_

> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_R2_cutadapt.txt` (26596)
<pre>
</pre>
Command completed. Elapsed time: 0:01:30. Running peak memory: 1.994GB.  
  PID: 26596;	Command: cutadapt;	Return code: 0;	Memory used: 0.29GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (26773,26774)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 1.994GB.  
  PID: 26773;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 26774;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	31605871	PEPPRO	_RES_

> `Trim_loss_rate`	67.66	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (26836)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_1_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_1_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_1_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_1_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 1.994GB.  
  PID: 26836;	Command: fastqc;	Return code: 0;	Memory used: 0.169GB

> `FastQC report r1`	fastqc/H9_PRO-seq_1_R2_trimmed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.histogram`  

> `fastq_pair -t 87956874 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq` (26949)
<pre>
Left paired: 30692275		Right paired: 30692275
Left single: 173668		Right single: 1557501
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:03:35. Running peak memory: 6.451GB.  
  PID: 26949;	Command: fastq_pair;	Return code: 0;	Memory used: 6.451GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_noadap.fastq.paired.fq -o H9_PRO-seq_1 -d /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt` (27553)
<pre>
</pre>
Command completed. Elapsed time: 0:01:48. Running peak memory: 6.451GB.  
  PID: 27553;	Command: flash;	Return code: 0;	Memory used: 0.087GB


### Plot adapter insertion distribution (02-18 11:29:35) elapsed: 906.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt -u 8` (27958)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.451GB.  
  PID: 27958;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:29:41) elapsed: 6.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/cutadapt/H9_PRO-seq_1.hist`

> `Degradation_ratio`	0.8896	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq` (27985)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 6.451GB.  
  PID: 27985;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq.paired.fq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq.paired.fq`  

> `fastq_pair -t 87956874 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (28002)
<pre>
Left paired: 30003416		Right paired: 30003416
Left single: 217927		Right single: 1602455
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:46. Running peak memory: 6.451GB.  
  PID: 28002;	Command: fastq_pair;	Return code: 0;	Memory used: 5.834GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/repaired.flag`  

> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq` (28482)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.451GB.  
  PID: 28482;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq` (28483)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 28483;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/repaired.flag` (28484)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 28484;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq.paired.fq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq`  

> `fastq_pair -t 87956874 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq` (28486)
<pre>
Left paired: 26674489		Right paired: 26674489
Left single: 126830		Right single: 4931382
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:18. Running peak memory: 6.451GB.  
  PID: 28486;	Command: fastq_pair;	Return code: 0;	Memory used: 6.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/dups_repaired.flag`  

> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq` (28835)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 28835;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq` (28837)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 28837;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/dups_repaired.flag` (28838)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 28838;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-18 11:35:05) elapsed: 324.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:35:05) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2` (28840)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 28840;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq` (28841)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
30003416 reads; of these:
  30003416 (100.00%) were unpaired; of these:
    27299076 (90.99%) aligned 0 times
    2704340 (9.01%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.01% overall alignment rate

> `Aligned_reads_human_rDNA`	5408680.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	17.11	PEPPRO	_RES_

### Map to human_rDNA (02-18 11:40:30) elapsed: 325.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2` (29447)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 29447;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R2.fq` (29449)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/fastq/H9_PRO-seq_1_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2704340 reads skipped
0 reads lost

### Map to genome (02-18 11:44:27) elapsed: 237.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_1 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/tmpnp63l3q8 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam` (29717,29718,29719)
<pre>
2176982 reads skipped
0 reads lost
27299076 reads; of these:
  27299076 (100.00%) were paired; of these:
    9566070 (35.04%) aligned concordantly 0 times
    14769056 (54.10%) aligned concordantly exactly 1 time
    2963950 (10.86%) aligned concordantly >1 times
    ----
    9566070 pairs aligned concordantly 0 times; of these:
      2129010 (22.26%) aligned discordantly 1 time
    ----
    7437060 pairs aligned 0 times concordantly or discordantly; of these:
      14874120 mates make up the pairs; of these:
        6258301 (42.08%) aligned 0 times
        3282618 (22.07%) aligned exactly 1 time
        5333201 (35.86%) aligned >1 times
88.54% overall alignment rate
[bam_sort_core] merging from 14 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 1:15:04. Running peak memory: 6.451GB.  
  PID: 29718;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 29717;	Command: bowtie2;	Return code: 0;	Memory used: 3.69GB  
  PID: 29719;	Command: samtools;	Return code: 0;	Memory used: 0.873GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (37858)
<pre>
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 6.451GB.  
  PID: 37858;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	48339851	PEPPRO	_RES_

> `QC_filtered_reads`	27585022	PEPPRO	_RES_

> `Aligned_reads`	20754829.0	PEPPRO	_RES_

> `Alignment_rate`	65.67	PEPPRO	_RES_

> `Total_efficiency`	21.24	PEPPRO	_RES_

> `Read_depth`	3.26	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_1 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/tmpnp63l3q8 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam` (40362,40363,40364)
<pre>
24497507 reads; of these:
  24497507 (100.00%) were paired; of these:
    8334941 (34.02%) aligned concordantly 0 times
    13482466 (55.04%) aligned concordantly exactly 1 time
    2680100 (10.94%) aligned concordantly >1 times
    ----
    8334941 pairs aligned concordantly 0 times; of these:
      1950931 (23.41%) aligned discordantly 1 time
    ----
    6384010 pairs aligned 0 times concordantly or discordantly; of these:
      12768020 mates make up the pairs; of these:
        5376104 (42.11%) aligned 0 times
        3000379 (23.50%) aligned exactly 1 time
        4391537 (34.39%) aligned >1 times
89.03% overall alignment rate
[bam_sort_core] merging from 13 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 1:09:19. Running peak memory: 6.451GB.  
  PID: 40363;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 40362;	Command: bowtie2;	Return code: 0;	Memory used: 3.665GB  
  PID: 40364;	Command: samtools;	Return code: 0;	Memory used: 0.93GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam` (49988)
<pre>
</pre>
Command completed. Elapsed time: 0:02:03. Running peak memory: 6.451GB.  
  PID: 49988;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (02-18 14:27:54) elapsed: 9807.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R1.fq` (50752)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 6.451GB.  
  PID: 50752;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/prealignments/H9_PRO-seq_1_human_rDNA_unmap_R2.fq` (50843)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 6.451GB.  
  PID: 50843;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam` (50915)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 6.451GB.  
  PID: 50915;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	674956	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (50996)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.451GB.  
  PID: 50996;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_noMT.bam` (51285,51286,51287,51288)
<pre>
</pre>
Command completed. Elapsed time: 0:00:43. Running peak memory: 6.451GB.  
  PID: 51287;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 51285;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 51286;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 51288;	Command: xargs;	Return code: 0;	Memory used: 0.083GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (51896)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 51896;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam` (51897)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.451GB.  
  PID: 51897;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


### Split BAM file (02-18 14:31:32) elapsed: 217.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam` (51934,51935)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:03:22. Running peak memory: 6.451GB.  
  PID: 51934;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 51935;	Command: samtools;	Return code: 0;	Memory used: 5.302GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE2.bam` (52604,52605)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:36. Running peak memory: 6.451GB.  
  PID: 52604;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 52605;	Command: samtools;	Return code: 0;	Memory used: 4.424GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam` (53613)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 6.451GB.  
  PID: 53613;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam` (53655,53656)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:54. Running peak memory: 6.451GB.  
  PID: 53655;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 53656;	Command: samtools;	Return code: 0;	Memory used: 4.909GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE2.bam` (54261,54262)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:37. Running peak memory: 6.451GB.  
  PID: 54261;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 54262;	Command: samtools;	Return code: 0;	Memory used: 4.089GB


### Calculate library complexity (02-18 14:44:22) elapsed: 770.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_out.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam` (54593)
<pre>
BAM_INPUT
TOTAL READS     = 20330230
COUNTS_SUM      = 20330230
DISTINCT READS  = 1.66542e+07
DISTINCT COUNTS = 262
MAX COUNT       = 16193
COUNTS OF 1     = 1.47277e+07
OBSERVED COUNTS (16194)
1	14727710
2	1349104
3	301309
4	112637
5	55354
6	31165
7	18956
8	12770
9	8912
10	6479
11	4906
12	3639
13	2910
14	2335
15	1909
16	1594
17	1295
18	1055
19	904
20	791
21	668
22	611
23	529
24	477
25	427
26	388
27	364
28	323
29	242
30	234
31	263
32	232
33	207
34	182
35	165
36	158
37	167
38	119
39	133
40	114
41	111
42	122
43	105
44	82
45	73
46	82
47	75
48	75
49	76
50	60
51	49
52	59
53	63
54	55
55	46
56	48
57	54
58	37
59	41
60	42
61	47
62	48
63	35
64	35
65	32
66	27
67	30
68	32
69	27
70	18
71	21
72	23
73	16
74	16
75	21
76	17
77	15
78	13
79	24
80	11
81	16
82	21
83	21
84	11
85	11
86	10
87	23
88	12
89	16
90	16
91	11
92	12
93	10
94	8
95	9
96	12
97	16
98	6
99	14
100	3
101	5
102	7
103	9
104	5
105	11
106	10
107	7
108	5
109	7
110	7
111	7
112	4
113	7
114	3
115	4
116	11
117	9
118	4
119	5
122	7
123	3
124	3
125	8
126	3
127	4
128	4
129	3
130	3
131	7
132	6
133	4
134	4
135	7
136	7
137	2
138	5
139	4
140	3
141	1
142	4
143	4
144	2
145	3
146	4
147	5
148	2
149	4
150	1
151	1
152	1
154	1
155	5
156	3
158	1
161	1
162	3
163	2
164	1
165	3
166	3
167	1
169	1
170	2
171	1
172	1
173	2
174	1
175	1
176	3
178	2
179	4
181	1
182	1
183	1
184	1
186	1
187	2
188	1
190	1
192	2
193	1
196	1
198	1
199	1
200	1
202	1
203	2
204	2
205	1
206	1
207	1
208	1
209	1
212	1
214	1
216	3
220	1
221	1
222	1
223	1
224	2
227	1
228	1
232	1
233	1
235	1
237	2
241	1
249	1
251	1
253	1
254	1
255	1
256	1
257	1
258	2
262	1
264	1
266	1
268	1
269	1
275	1
277	1
280	1
296	1
299	1
303	1
306	1
313	1
315	1
317	2
325	1
345	1
348	1
406	1
412	1
430	1
433	1
436	1
437	1
493	1
494	1
540	1
547	1
606	2
611	1
612	1
734	1
787	1
789	1
848	1
927	1
943	1
1376	1
1390	1
1414	1
2008	1
4768	1
5683	1
6998	1
14724	1
16193	1

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
sample size: 20000000
</pre>
Command completed. Elapsed time: 0:01:49. Running peak memory: 6.451GB.  
  PID: 54593;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam` (54917)
<pre>
BAM_INPUT
TOTAL READS     = 20330230
DISTINCT READS  = 1.66542e+07
DISTINCT COUNTS = 262
MAX COUNT       = 16193
COUNTS OF 1     = 1.47277e+07
MAX TERMS       = 100
OBSERVED COUNTS (16194)
1	14727710
2	1349104
3	301309
4	112637
5	55354
6	31165
7	18956
8	12770
9	8912
10	6479
11	4906
12	3639
13	2910
14	2335
15	1909
16	1594
17	1295
18	1055
19	904
20	791
21	668
22	611
23	529
24	477
25	427
26	388
27	364
28	323
29	242
30	234
31	263
32	232
33	207
34	182
35	165
36	158
37	167
38	119
39	133
40	114
41	111
42	122
43	105
44	82
45	73
46	82
47	75
48	75
49	76
50	60
51	49
52	59
53	63
54	55
55	46
56	48
57	54
58	37
59	41
60	42
61	47
62	48
63	35
64	35
65	32
66	27
67	30
68	32
69	27
70	18
71	21
72	23
73	16
74	16
75	21
76	17
77	15
78	13
79	24
80	11
81	16
82	21
83	21
84	11
85	11
86	10
87	23
88	12
89	16
90	16
91	11
92	12
93	10
94	8
95	9
96	12
97	16
98	6
99	14
100	3
101	5
102	7
103	9
104	5
105	11
106	10
107	7
108	5
109	7
110	7
111	7
112	4
113	7
114	3
115	4
116	11
117	9
118	4
119	5
122	7
123	3
124	3
125	8
126	3
127	4
128	4
129	3
130	3
131	7
132	6
133	4
134	4
135	7
136	7
137	2
138	5
139	4
140	3
141	1
142	4
143	4
144	2
145	3
146	4
147	5
148	2
149	4
150	1
151	1
152	1
154	1
155	5
156	3
158	1
161	1
162	3
163	2
164	1
165	3
166	3
167	1
169	1
170	2
171	1
172	1
173	2
174	1
175	1
176	3
178	2
179	4
181	1
182	1
183	1
184	1
186	1
187	2
188	1
190	1
192	2
193	1
196	1
198	1
199	1
200	1
202	1
203	2
204	2
205	1
206	1
207	1
208	1
209	1
212	1
214	1
216	3
220	1
221	1
222	1
223	1
224	2
227	1
228	1
232	1
233	1
235	1
237	2
241	1
249	1
251	1
253	1
254	1
255	1
256	1
257	1
258	2
262	1
264	1
266	1
268	1
269	1
275	1
277	1
280	1
296	1
299	1
303	1
306	1
313	1
315	1
317	2
325	1
345	1
348	1
406	1
412	1
430	1
433	1
436	1
437	1
493	1
494	1
540	1
547	1
606	2
611	1
612	1
734	1
787	1
789	1
848	1
927	1
943	1
1376	1
1390	1
1414	1
2008	1
4768	1
5683	1
6998	1
14724	1
16193	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
......................._.............................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:58. Running peak memory: 6.451GB.  
  PID: 54917;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_counts.txt` (55191)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.451GB.  
  PID: 55191;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_plot` (55248)
<pre>
Processing H9_PRO-seq_1
INFO: Found real counts for H9_PRO-seq_1 - Total (M): 21.966708 Unique (M): 20.33023

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.451GB.  
  PID: 55248;	Command: Rscript;	Return code: 0;	Memory used: 0.237GB

> `Library complexity`	QC_hg38/H9_PRO-seq_1_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_1_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8774	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 14:48:45) elapsed: 263.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam` (55267)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 6.451GB.  
  PID: 55267;	Command: samtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv` (55282)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/tmp_H9_PRO-seq_1_PE1_1jkyqzj3'
Processing with 8 cores...
Discarding 96 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270720v1_random', 'chr14_KI270724v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 99 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270539v1', 'chrUn_KI270587v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.451GB.  
  PID: 55282;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.236GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	10983354.0	PEPPRO	_RES_

> `PBC2`	10983354.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_unmap.bam` (55336)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 6.451GB.  
  PID: 55336;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_temp.bam`

> `Unmapped_reads`	6258301	PEPPRO	_RES_

### Split BAM by strand (02-18 14:49:51) elapsed: 66.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam` (55381)
<pre>
</pre>
Command completed. Elapsed time: 0:01:17. Running peak memory: 6.451GB.  
  PID: 55381;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam` (55696)
<pre>
</pre>
Command completed. Elapsed time: 0:01:14. Running peak memory: 6.451GB.  
  PID: 55696;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 14:52:22) elapsed: 151.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (55919)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 55919;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt` (55921)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.451GB.  
  PID: 55921;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.513GB


> `TSS_Plus_Score`	26.7	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt` (55952)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.451GB.  
  PID: 55952;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.508GB


> `TSS_Minus_Score`	10.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt` (55980)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.451GB.  
  PID: 55980;	Command: Rscript;	Return code: 0;	Memory used: 0.239GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_1_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_1_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt` (55999,56000,56001,56002)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 55999;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 56001;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 56000;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 56002;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt` (56004)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 56004;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-18 14:52:44) elapsed: 22.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed` (56006,56007)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.451GB.  
  PID: 56006;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 56007;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed` (56011,56012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 56011;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 56012;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSS_density.bed` (56014,56015,56016,56017)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 6.451GB.  
  PID: 56015;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 56017;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 56014;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 56016;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_gene_body_density.bed` (56047,56048,56049)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 6.451GB.  
  PID: 56047;	Command: bedtools;	Return code: 0;	Memory used: 0.09GB  
  PID: 56049;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 56048;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed` (56094,56095,56096)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 56094;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 56096;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 56095;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	32.63	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed` (56101)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.451GB.  
  PID: 56101;	Command: Rscript;	Return code: 0;	Memory used: 0.238GB

> `Pause index`	QC_hg38/H9_PRO-seq_1_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_1_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_pause_index.bed` (56123)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 56123;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 14:54:06) elapsed: 82.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam`
20754829.0 7959339

> `Plus_FRiP`	0.38	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam`
20754829.0 7448642

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_gene_sort.bed` (56185,56186)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.451GB.  
  PID: 56185;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 56186;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_gene_coverage.bed` (56189)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 6.451GB.  
  PID: 56189;	Command: bedtools;	Return code: 0;	Memory used: 0.07GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed.gz` (56437)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 56437;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed` (56439)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 56439;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 14:55:33) elapsed: 88.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/raw/hg38_annotations.bed` (56447)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.451GB.  
  PID: 56447;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR"` (56449)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 56449;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed` (56450,56451,56452,56453)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.451GB.  
  PID: 56450;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 56451;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 56453;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 56452;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed` (56456)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.451GB.  
  PID: 56456;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_minus_coverage.bed` (56515)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.451GB.  
  PID: 56515;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR"` (56531)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 56531;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed` (56532,56533,56534,56535)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.451GB.  
  PID: 56532;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 56533;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 56535;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 56534;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed` (56537)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.451GB.  
  PID: 56537;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_minus_coverage.bed` (56581)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.451GB.  
  PID: 56581;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed` (56677,56678,56679,56680)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.451GB.  
  PID: 56677;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 56678;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 56680;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB  
  PID: 56679;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed` (56682)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.451GB.  
  PID: 56682;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_minus_coverage.bed` (56712)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.451GB.  
  PID: 56712;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed` (56725,56726,56727,56728)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.451GB.  
  PID: 56725;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 56726;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 56728;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 56727;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed` (56735)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.451GB.  
  PID: 56735;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_minus_coverage.bed` (56750)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.451GB.  
  PID: 56750;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed` (56773,56774,56775,56776)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.451GB.  
  PID: 56773;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 56775;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 56774;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 56776;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed` (56779)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 6.451GB.  
  PID: 56779;	Command: bedtools;	Return code: 0;	Memory used: 0.055GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_minus_coverage.bed` (56797)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 6.451GB.  
  PID: 56797;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed` (56817,56818,56819,56820)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 56817;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 56818;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 56820;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 56819;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed` (56822)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 6.451GB.  
  PID: 56822;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_minus_coverage.bed` (56846)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.451GB.  
  PID: 56846;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region"` (56861)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 56861;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed` (56862,56863,56864,56865)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.451GB.  
  PID: 56862;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 56864;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 56863;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 56865;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (56868)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.451GB.  
  PID: 56868;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed` (56884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.451GB.  
  PID: 56884;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB


### Plot FRiF/PRiF (02-18 14:59:20) elapsed: 226.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_1 -z 3099922541 -n 11202246 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (56910)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 6.451GB.  
  PID: 56910;	Command: Rscript;	Return code: 0;	Memory used: 0.463GB

> `FRiF`	QC_hg38/H9_PRO-seq_1_frif.pdf	FRiF	QC_hg38/H9_PRO-seq_1_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_1 -z 3099922541 -n 11202246 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (56959)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.451GB.  
  PID: 56959;	Command: Rscript;	Return code: 0;	Memory used: 0.463GB

> `PRiF`	QC_hg38/H9_PRO-seq_1_prif.pdf	PRiF	QC_hg38/H9_PRO-seq_1_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 15:00:23) elapsed: 64.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_exons_sort.bed` (57229,57230)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.451GB.  
  PID: 57230;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB  
  PID: 57229;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_introns_sort.bed` (57241,57242,57243)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.451GB.  
  PID: 57241;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 57243;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 57242;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_coverage.bed` (57266)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.451GB.  
  PID: 57266;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_coverage.bed` (57335)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 6.451GB.  
  PID: 57335;	Command: bedtools;	Return code: 0;	Memory used: 0.044GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/20.754829)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_rpkm.bed` (57495,57496,57497)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.451GB.  
  PID: 57495;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 57497;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 57496;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/20.754829)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_rpkm.bed` (57500,57501,57502)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.451GB.  
  PID: 57500;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 57502;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 57501;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed` (57504,57505,57506)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 57504;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 57506;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 57505;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.23	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed --annotate` (57512)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.451GB.  
  PID: 57512;	Command: Rscript;	Return code: 0;	Memory used: 0.321GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_1_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_1_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/QC_hg38/H9_PRO-seq_1_exon_intron_ratios.bed` (57531)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.451GB.  
  PID: 57531;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-18 15:01:42) elapsed: 79.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam` (57539)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.451GB.  
  PID: 57539;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (57580)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_1_plus_cuttrace_ai7s1lg3'
Processing with 2 cores...
stdin is empty of data
Discarding 115 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr14_KI270724v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 80 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270438v1', 'chrUn_KI270515v1', 'chrUn_KI270539v1', 'chrUn_KI270587v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 80 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_exact_body_0-mer.bw'
Merging 80 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:34. Running peak memory: 6.451GB.  
  PID: 57580;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.139GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam` (59656)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.451GB.  
  PID: 59656;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_exact_body_0-mer.bw -p 5 --variable-step --tail-edge` (59670)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/aligned_hg38/H9_PRO-seq_1_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_1_minus_cuttrace_kezy_tzk'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_KI270722v1_random', 'chr14_KI270724v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_1/signal_hg38/H9_PRO-seq_1_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 6.451GB.  
  PID: 59670;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.272GB

Starting cleanup: 68 files; 4 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  4:02:57
*  Total elapsed time (all runs):  8:16:53
*         Peak memory (this run):  6.4513 GB
*        Pipeline completed time: 2020-02-18 15:11:55
