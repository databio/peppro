### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name H9_PRO-seq_3 --genome hg38 --input /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz --protocol PRO --umi-len 8 --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj40-15c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/
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
*              `input`:  `['/project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_3`
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

Local input file: /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz
Local input file: /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz

> `File_mb`	2340.86	PEPPRO	_RES_

> `Read_type`	paired	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:08:59) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE1.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz` (53172)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 53172;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep3_PE2.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz` (53176)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 53176;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq` (53180)
<pre>
</pre>
Command completed. Elapsed time: 0:02:39. Running peak memory: 0.003GB.  
  PID: 53180;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2.fastq` (271823)
<pre>
</pre>
Command completed. Elapsed time: 0:01:40. Running peak memory: 0.003GB.  
  PID: 271823;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	92628962	PEPPRO	_RES_

> `Fastq_reads`	92628962	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/H9_PRO-seq_3_R2.fastq.gz']

### FASTQ processing:  (02-18 11:14:13) elapsed: 314.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt` (107756)
<pre>
</pre>
Command completed. Elapsed time: 0:01:21. Running peak memory: 0.272GB.  
  PID: 107756;	Command: cutadapt;	Return code: 0;	Memory used: 0.272GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (25636,25637)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 0.272GB.  
  PID: 25636;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 25637;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	27419170	PEPPRO	_RES_

> `Trim_loss_rate`	70.4	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (25714)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_3_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_3_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_3_R1_processed.fastq
Analysis complete for H9_PRO-seq_3_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 0.272GB.  
  PID: 25714;	Command: fastqc;	Return code: 0;	Memory used: 0.177GB

> `FastQC report r1`	fastqc/H9_PRO-seq_3_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_dedup.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq` (25967)
<pre>
[INFO][0m 1862033 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 1.993GB.  
  PID: 25967;	Command: seqkit;	Return code: 0;	Memory used: 1.993GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq` (26060,26061)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 1.993GB.  
  PID: 26060;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 26061;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	34725419.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	18106301.0	PEPPRO	_RES_

> `Duplicate_reads`	1862033.0	PEPPRO	_RES_

> `Pct_reads_too_short`	39.0943	PEPPRO	_RES_

> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_R2_cutadapt.txt` (26168)
<pre>
</pre>
Command completed. Elapsed time: 0:01:22. Running peak memory: 1.993GB.  
  PID: 26168;	Command: cutadapt;	Return code: 0;	Memory used: 0.288GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (26530,26531)
<pre>
</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 1.993GB.  
  PID: 26530;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 26531;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	28703909	PEPPRO	_RES_

> `Trim_loss_rate`	69.01	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (26633)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_3_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_3_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_3_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_3_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 1.993GB.  
  PID: 26633;	Command: fastqc;	Return code: 0;	Memory used: 0.172GB

> `FastQC report r1`	fastqc/H9_PRO-seq_3_R2_trimmed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.histogram`  

> `fastq_pair -t 83366065 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq` (26764)
<pre>
Left paired: 28081513		Right paired: 28081513
Left single: 126667		Right single: 1382160
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:48. Running peak memory: 5.443GB.  
  PID: 26764;	Command: fastq_pair;	Return code: 0;	Memory used: 5.443GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_noadap.fastq.paired.fq -o H9_PRO-seq_3 -d /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt` (27003)
<pre>
</pre>
Command completed. Elapsed time: 0:01:23. Running peak memory: 5.443GB.  
  PID: 27003;	Command: flash;	Return code: 0;	Memory used: 0.089GB


### Plot adapter insertion distribution (02-18 11:26:23) elapsed: 730.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt -u 8` (27434)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 5.443GB.  
  PID: 27434;	Command: Rscript;	Return code: 0;	Memory used: 0.368GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_3_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_3_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:26:37) elapsed: 14.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/cutadapt/H9_PRO-seq_3.hist`

> `Degradation_ratio`	0.9769	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq` (27476)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 5.443GB.  
  PID: 27476;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq.paired.fq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq.paired.fq`  

> `fastq_pair -t 83366065 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (27500)
<pre>
Left paired: 27272197		Right paired: 27272197
Left single: 146973		Right single: 1431712
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:10. Running peak memory: 5.688GB.  
  PID: 27500;	Command: fastq_pair;	Return code: 0;	Memory used: 5.688GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/repaired.flag`  

> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq` (27921)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.688GB.  
  PID: 27921;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq` (27922)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 27922;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/repaired.flag` (27924)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 27924;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq.paired.fq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq.paired.fq`  

> `fastq_pair -t 83366065 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq` (27925)
<pre>
Left paired: 25551175		Right paired: 25551175
Left single: 94450		Right single: 3152734
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:19. Running peak memory: 5.688GB.  
  PID: 27925;	Command: fastq_pair;	Return code: 0;	Memory used: 5.139GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/dups_repaired.flag`  

> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq` (28353)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 28353;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq` (28355)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.688GB.  
  PID: 28355;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/dups_repaired.flag` (28357)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 28357;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-18 11:31:30) elapsed: 292.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:31:30) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2` (28358)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 28358;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq` (28359)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_3 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
27272197 reads; of these:
  27272197 (100.00%) were unpaired; of these:
    24569260 (90.09%) aligned 0 times
    2702937 (9.91%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.91% overall alignment rate

> `Aligned_reads_human_rDNA`	5405874.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	18.83	PEPPRO	_RES_

### Map to human_rDNA (02-18 11:34:34) elapsed: 184.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2` (28605)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 28605;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R2.fq` (28606)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_3 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/fastq/H9_PRO-seq_3_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2702937 reads skipped
0 reads lost

### Map to genome (02-18 11:37:34) elapsed: 181.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_3 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/tmp3p3e5792 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam` (29032,29033,29034)
<pre>
2244415 reads skipped
0 reads lost
24569260 reads; of these:
  24569260 (100.00%) were paired; of these:
    9705626 (39.50%) aligned concordantly 0 times
    12263897 (49.92%) aligned concordantly exactly 1 time
    2599737 (10.58%) aligned concordantly >1 times
    ----
    9705626 pairs aligned concordantly 0 times; of these:
      2262313 (23.31%) aligned discordantly 1 time
    ----
    7443313 pairs aligned 0 times concordantly or discordantly; of these:
      14886626 mates make up the pairs; of these:
        5843431 (39.25%) aligned 0 times
        3306807 (22.21%) aligned exactly 1 time
        5736388 (38.53%) aligned >1 times
88.11% overall alignment rate
[bam_sort_core] merging from 13 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 1:01:26. Running peak memory: 5.688GB.  
  PID: 29033;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 29032;	Command: bowtie2;	Return code: 0;	Memory used: 3.661GB  
  PID: 29034;	Command: samtools;	Return code: 0;	Memory used: 0.874GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (35657)
<pre>
</pre>
Command completed. Elapsed time: 0:01:27. Running peak memory: 5.688GB.  
  PID: 35657;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	43295089	PEPPRO	_RES_

> `QC_filtered_reads`	25127415	PEPPRO	_RES_

> `Aligned_reads`	18167674.5	PEPPRO	_RES_

> `Alignment_rate`	63.29	PEPPRO	_RES_

> `Total_efficiency`	19.61	PEPPRO	_RES_

> `Read_depth`	2.98	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_3 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/tmp3p3e5792 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam` (37212,37217,37218)
<pre>
23306760 reads; of these:
  23306760 (100.00%) were paired; of these:
    9052177 (38.84%) aligned concordantly 0 times
    11774634 (50.52%) aligned concordantly exactly 1 time
    2479949 (10.64%) aligned concordantly >1 times
    ----
    9052177 pairs aligned concordantly 0 times; of these:
      2177215 (24.05%) aligned discordantly 1 time
    ----
    6874962 pairs aligned 0 times concordantly or discordantly; of these:
      13749924 mates make up the pairs; of these:
        5423690 (39.45%) aligned 0 times
        3177751 (23.11%) aligned exactly 1 time
        5148483 (37.44%) aligned >1 times
88.36% overall alignment rate
[bam_sort_core] merging from 12 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:56:57. Running peak memory: 5.688GB.  
  PID: 37212;	Command: bowtie2;	Return code: 0;	Memory used: 3.66GB  
  PID: 37217;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 37218;	Command: samtools;	Return code: 0;	Memory used: 0.873GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam` (45041)
<pre>
</pre>
Command completed. Elapsed time: 0:01:21. Running peak memory: 5.688GB.  
  PID: 45041;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (02-18 13:51:55) elapsed: 8061.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R2.fq` (45137)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 5.688GB.  
  PID: 45137;	Command: pigz;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/prealignments/H9_PRO-seq_3_human_rDNA_unmap_R1.fq` (45177)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 5.688GB.  
  PID: 45177;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam` (45219)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 5.688GB.  
  PID: 45219;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	711009	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (45265)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 5.688GB.  
  PID: 45265;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_noMT.bam` (45288,45289,45290,45291)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 5.688GB.  
  PID: 45290;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 45288;	Command: samtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 45289;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 45291;	Command: xargs;	Return code: 0;	Memory used: 0.061GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (45360)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 45360;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam` (45361)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 5.688GB.  
  PID: 45361;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


### Split BAM file (02-18 13:55:10) elapsed: 195.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam` (45594,45595)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:42. Running peak memory: 5.688GB.  
  PID: 45594;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 45595;	Command: samtools;	Return code: 0;	Memory used: 4.686GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE2.bam` (45819,45820)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:15. Running peak memory: 5.688GB.  
  PID: 45819;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 45820;	Command: samtools;	Return code: 0;	Memory used: 3.785GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam` (46295)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 5.688GB.  
  PID: 46295;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam` (46343,46344)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:38. Running peak memory: 5.688GB.  
  PID: 46343;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 46344;	Command: samtools;	Return code: 0;	Memory used: 4.574GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE2.bam` (46616,46617)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:10. Running peak memory: 5.688GB.  
  PID: 46616;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 46617;	Command: samtools;	Return code: 0;	Memory used: 3.693GB


### Calculate library complexity (02-18 14:06:11) elapsed: 661.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_out.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam` (46991)
<pre>
BAM_INPUT
TOTAL READS     = 19005606
COUNTS_SUM      = 19005606
DISTINCT READS  = 1.58447e+07
DISTINCT COUNTS = 285
MAX COUNT       = 19581
COUNTS OF 1     = 1.42036e+07
OBSERVED COUNTS (19582)
1	14203580
2	1156600
3	254232
4	94789
5	45836
6	26066
7	16009
8	10530
9	7388
10	5314
11	3937
12	3035
13	2237
14	1796
15	1517
16	1252
17	1031
18	892
19	767
20	636
21	535
22	462
23	411
24	386
25	367
26	297
27	272
28	227
29	236
30	215
31	209
32	162
33	160
34	160
35	150
36	126
37	130
38	108
39	116
40	107
41	94
42	95
43	81
44	70
45	69
46	67
47	68
48	59
49	70
50	59
51	60
52	65
53	42
54	43
55	50
56	44
57	45
58	37
59	39
60	31
61	47
62	32
63	28
64	33
65	38
66	28
67	31
68	33
69	29
70	31
71	21
72	20
73	18
74	21
75	24
76	18
77	28
78	20
79	17
80	17
81	14
82	17
83	14
84	13
85	13
86	17
87	11
88	17
89	14
90	8
91	13
92	10
93	15
94	14
95	11
96	5
97	13
98	10
99	10
100	12
101	3
102	14
103	11
104	8
105	15
106	9
107	9
108	9
109	6
110	9
111	5
112	4
113	4
114	6
115	12
116	1
117	7
118	6
119	9
120	4
121	2
122	3
123	4
124	5
125	4
126	3
127	4
128	7
129	7
130	2
131	5
132	3
133	5
134	8
135	3
136	5
137	4
138	7
139	2
140	2
141	3
142	4
143	3
144	6
145	3
146	6
147	2
148	1
149	5
150	2
151	3
152	4
153	5
154	2
155	5
156	3
157	4
158	5
159	2
160	4
161	2
162	1
163	3
164	1
166	3
167	1
168	3
169	3
170	3
172	3
173	6
174	1
175	2
176	2
177	1
178	2
179	1
181	1
183	1
184	3
185	3
186	1
187	1
188	2
189	1
190	1
191	2
193	3
194	1
195	2
196	2
197	3
199	3
200	2
201	2
202	2
203	4
204	1
205	1
207	4
208	1
210	1
211	1
212	1
213	2
214	2
216	1
217	2
223	3
224	1
225	1
226	1
229	2
231	2
233	1
234	1
235	1
238	1
239	1
242	2
246	1
247	1
248	3
250	1
254	1
256	2
259	2
260	1
264	1
271	1
274	1
276	1
284	2
294	1
299	1
304	1
305	1
318	1
334	1
335	1
348	1
351	1
353	1
357	1
366	1
390	2
391	1
393	1
394	1
397	1
422	1
432	1
443	1
450	1
467	1
469	1
500	1
508	1
513	1
516	1
526	1
537	1
543	1
560	1
577	1
582	1
611	1
619	1
653	1
697	1
887	1
939	1
964	1
1000	1
1280	1
1283	1
1718	1
1790	1
1929	1
2740	1
6648	1
7571	1
9539	1
16742	1
19581	1

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
</pre>
Command completed. Elapsed time: 0:01:44. Running peak memory: 5.688GB.  
  PID: 46991;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam` (47124)
<pre>
BAM_INPUT
TOTAL READS     = 19005606
DISTINCT READS  = 1.58447e+07
DISTINCT COUNTS = 285
MAX COUNT       = 19581
COUNTS OF 1     = 1.42036e+07
MAX TERMS       = 100
OBSERVED COUNTS (19582)
1	14203580
2	1156600
3	254232
4	94789
5	45836
6	26066
7	16009
8	10530
9	7388
10	5314
11	3937
12	3035
13	2237
14	1796
15	1517
16	1252
17	1031
18	892
19	767
20	636
21	535
22	462
23	411
24	386
25	367
26	297
27	272
28	227
29	236
30	215
31	209
32	162
33	160
34	160
35	150
36	126
37	130
38	108
39	116
40	107
41	94
42	95
43	81
44	70
45	69
46	67
47	68
48	59
49	70
50	59
51	60
52	65
53	42
54	43
55	50
56	44
57	45
58	37
59	39
60	31
61	47
62	32
63	28
64	33
65	38
66	28
67	31
68	33
69	29
70	31
71	21
72	20
73	18
74	21
75	24
76	18
77	28
78	20
79	17
80	17
81	14
82	17
83	14
84	13
85	13
86	17
87	11
88	17
89	14
90	8
91	13
92	10
93	15
94	14
95	11
96	5
97	13
98	10
99	10
100	12
101	3
102	14
103	11
104	8
105	15
106	9
107	9
108	9
109	6
110	9
111	5
112	4
113	4
114	6
115	12
116	1
117	7
118	6
119	9
120	4
121	2
122	3
123	4
124	5
125	4
126	3
127	4
128	7
129	7
130	2
131	5
132	3
133	5
134	8
135	3
136	5
137	4
138	7
139	2
140	2
141	3
142	4
143	3
144	6
145	3
146	6
147	2
148	1
149	5
150	2
151	3
152	4
153	5
154	2
155	5
156	3
157	4
158	5
159	2
160	4
161	2
162	1
163	3
164	1
166	3
167	1
168	3
169	3
170	3
172	3
173	6
174	1
175	2
176	2
177	1
178	2
179	1
181	1
183	1
184	3
185	3
186	1
187	1
188	2
189	1
190	1
191	2
193	3
194	1
195	2
196	2
197	3
199	3
200	2
201	2
202	2
203	4
204	1
205	1
207	4
208	1
210	1
211	1
212	1
213	2
214	2
216	1
217	2
223	3
224	1
225	1
226	1
229	2
231	2
233	1
234	1
235	1
238	1
239	1
242	2
246	1
247	1
248	3
250	1
254	1
256	2
259	2
260	1
264	1
271	1
274	1
276	1
284	2
294	1
299	1
304	1
305	1
318	1
334	1
335	1
348	1
351	1
353	1
357	1
366	1
390	2
391	1
393	1
394	1
397	1
422	1
432	1
443	1
450	1
467	1
469	1
500	1
508	1
513	1
516	1
526	1
537	1
543	1
560	1
577	1
582	1
611	1
619	1
653	1
697	1
887	1
939	1
964	1
1000	1
1280	1
1283	1
1718	1
1790	1
1929	1
2740	1
6648	1
7571	1
9539	1
16742	1
19581	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
.............................................._......................._.................._.............
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:51. Running peak memory: 5.688GB.  
  PID: 47124;	Command: preseq;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_counts.txt` (47257)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 5.688GB.  
  PID: 47257;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_plot` (47507)
<pre>
Processing H9_PRO-seq_3
INFO: Found real counts for H9_PRO-seq_3 - Total (M): 19.481041 Unique (M): 19.005606

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 5.688GB.  
  PID: 47507;	Command: Rscript;	Return code: 0;	Memory used: 0.292GB

> `Library complexity`	QC_hg38/H9_PRO-seq_3_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_3_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8818	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 14:10:19) elapsed: 248.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam` (47526)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 5.688GB.  
  PID: 47526;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv` (47540)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/tmp_H9_PRO-seq_3_PE1_fgx_iyuu'
Processing with 8 cores...
Discarding 101 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270748v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 94 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 5.688GB.  
  PID: 47540;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.123GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	9740520.5	PEPPRO	_RES_

> `PBC2`	9740520.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_unmap.bam` (47590)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 5.688GB.  
  PID: 47590;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_temp.bam`

> `Unmapped_reads`	5843431	PEPPRO	_RES_

### Split BAM by strand (02-18 14:11:10) elapsed: 52.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam` (47628)
<pre>
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 5.688GB.  
  PID: 47628;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam` (47700)
<pre>
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 5.688GB.  
  PID: 47700;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 14:13:28) elapsed: 138.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (47778)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 47778;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt` (47779)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 5.688GB.  
  PID: 47779;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.491GB


> `TSS_Plus_Score`	28.6	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt` (47809)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 5.688GB.  
  PID: 47809;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.487GB


> `TSS_Minus_Score`	10.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt` (47834)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 5.688GB.  
  PID: 47834;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_3_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_3_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt` (47863,47864,47865,47866)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 47863;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 47865;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 47864;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 47866;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt` (47868)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 47868;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (02-18 14:13:47) elapsed: 19.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed` (47870,47871)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 5.688GB.  
  PID: 47870;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 47871;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed` (47875,47876)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 47875;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 47876;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSS_density.bed` (47878,47879,47880,47881)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 5.688GB.  
  PID: 47879;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 47881;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 47878;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 47880;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_gene_body_density.bed` (47933,47934,47935)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 5.688GB.  
  PID: 47935;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 47933;	Command: bedtools;	Return code: 0;	Memory used: 0.062GB  
  PID: 47934;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed` (47999,48000,48001)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 47999;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 48001;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 48000;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	33.89	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed` (48007)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.688GB.  
  PID: 48007;	Command: Rscript;	Return code: 0;	Memory used: 0.238GB

> `Pause index`	QC_hg38/H9_PRO-seq_3_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_3_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_pause_index.bed` (48026)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 48026;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 14:14:56) elapsed: 68.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam`
18167674.5 7066591

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam`
18167674.5 6601482

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_gene_sort.bed` (48278,48279)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.688GB.  
  PID: 48278;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 48279;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_gene_coverage.bed` (48281)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 5.688GB.  
  PID: 48281;	Command: bedtools;	Return code: 0;	Memory used: 0.067GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed.gz` (48313)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 48313;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed` (48314)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 48314;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 14:16:02) elapsed: 67.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/raw/hg38_annotations.bed` (48323)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 5.688GB.  
  PID: 48323;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR"` (48325)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 48325;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed` (48326,48327,48328,48329)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.688GB.  
  PID: 48326;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 48327;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 48329;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 48328;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed` (48334)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.688GB.  
  PID: 48334;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_minus_coverage.bed` (48347)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 5.688GB.  
  PID: 48347;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR"` (48365)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 48365;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed` (48366,48367,48368,48369)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.688GB.  
  PID: 48366;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 48367;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 48369;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 48368;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed` (48371)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.688GB.  
  PID: 48371;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_minus_coverage.bed` (48387)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.688GB.  
  PID: 48387;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed` (48437,48438,48439,48440)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.688GB.  
  PID: 48437;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 48438;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 48440;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 48439;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed` (48442)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.688GB.  
  PID: 48442;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_minus_coverage.bed` (48477)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 5.688GB.  
  PID: 48477;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed` (48522,48523,48524,48525)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 5.688GB.  
  PID: 48522;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 48523;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 48525;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 48524;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed` (48532)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 5.688GB.  
  PID: 48532;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_minus_coverage.bed` (48551)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 5.688GB.  
  PID: 48551;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed` (48568,48569,48570,48571)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.688GB.  
  PID: 48568;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 48570;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 48569;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 48571;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed` (48574)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 5.688GB.  
  PID: 48574;	Command: bedtools;	Return code: 0;	Memory used: 0.06GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_minus_coverage.bed` (48593)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 5.688GB.  
  PID: 48593;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed` (48647,48648,48649,48650)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 48647;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 48648;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 48650;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 48649;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed` (48652)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.688GB.  
  PID: 48652;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_minus_coverage.bed` (48671)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.688GB.  
  PID: 48671;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region"` (48689)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 48689;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed` (48690,48692,48693,48694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.688GB.  
  PID: 48690;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 48693;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 48692;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 48694;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed` (48696)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.688GB.  
  PID: 48696;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_minus_coverage.bed` (48714)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 5.688GB.  
  PID: 48714;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


### Plot FRiF/PRiF (02-18 14:19:21) elapsed: 199.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_3 -z 3099922541 -n 9939223 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed` (48741)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 5.688GB.  
  PID: 48741;	Command: Rscript;	Return code: 0;	Memory used: 0.446GB

> `FRiF`	QC_hg38/H9_PRO-seq_3_frif.pdf	FRiF	QC_hg38/H9_PRO-seq_3_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_3 -z 3099922541 -n 9939223 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed` (48811)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 5.688GB.  
  PID: 48811;	Command: Rscript;	Return code: 0;	Memory used: 0.467GB

> `PRiF`	QC_hg38/H9_PRO-seq_3_prif.pdf	PRiF	QC_hg38/H9_PRO-seq_3_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 14:20:25) elapsed: 64.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_exons_sort.bed` (49081,49082)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.688GB.  
  PID: 49082;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB  
  PID: 49081;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_introns_sort.bed` (49102,49103,49104)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 5.688GB.  
  PID: 49102;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 49104;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 49103;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_coverage.bed` (49112)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 5.688GB.  
  PID: 49112;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_coverage.bed` (49144)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 5.688GB.  
  PID: 49144;	Command: bedtools;	Return code: 0;	Memory used: 0.058GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.1676745)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_rpkm.bed` (49187,49188,49189)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.688GB.  
  PID: 49187;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 49189;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 49188;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.1676745)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_rpkm.bed` (49192,49193,49194)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 5.688GB.  
  PID: 49192;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 49194;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 49193;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed` (49196,49197,49198)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 49196;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 49198;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 49197;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.28	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed --annotate` (49205)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 5.688GB.  
  PID: 49205;	Command: Rscript;	Return code: 0;	Memory used: 0.366GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_3_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_3_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/QC_hg38/H9_PRO-seq_3_exon_intron_ratios.bed` (49225)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 5.688GB.  
  PID: 49225;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-18 14:21:35) elapsed: 70.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam` (49233)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 5.688GB.  
  PID: 49233;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (49242)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_3_plus_cuttrace_04cwhflw'
Processing with 2 cores...
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:57. Running peak memory: 5.688GB.  
  PID: 49242;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.546GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam` (51956)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 5.688GB.  
  PID: 51956;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_exact_body_0-mer.bw -p 5 --variable-step --tail-edge` (51980)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/aligned_hg38/H9_PRO-seq_3_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_3_minus_cuttrace_k52mwomk'
Processing with 5 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270737v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_GL000226v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_3/signal_hg38/H9_PRO-seq_3_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 5.688GB.  
  PID: 51980;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.335GB

Starting cleanup: 68 files; 4 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:23:06
*  Total elapsed time (all runs):  6:50:53
*         Peak memory (this run):  5.6883 GB
*        Pipeline completed time: 2020-02-18 14:32:04
