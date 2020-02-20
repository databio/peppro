### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name H9_PRO-seq_2 --genome hg38 --input /project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz --protocol PRO --umi-len 8 --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj40-15c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/
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
*              `input`:  `['/project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_2`
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

Local input file: /project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz
Local input file: /project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz

> `File_mb`	2754.18	PEPPRO	_RES_

> `Read_type`	paired	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:08:59) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz` (53171)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 53171;	Command: ln;	Return code: 0;	Memory used: 0.002GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz` (53175)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 53175;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq` (53179)
<pre>
</pre>
Command completed. Elapsed time: 0:02:55. Running peak memory: 0.003GB.  
  PID: 53179;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2.fastq` (398411)
<pre>
</pre>
Command completed. Elapsed time: 0:01:55. Running peak memory: 0.003GB.  
  PID: 398411;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	115714026	PEPPRO	_RES_

> `Fastq_reads`	115714026	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz']

### FASTQ processing:  (02-18 11:15:54) elapsed: 415.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt` (25671)
<pre>
</pre>
Command completed. Elapsed time: 0:01:47. Running peak memory: 0.279GB.  
  PID: 25671;	Command: cutadapt;	Return code: 0;	Memory used: 0.279GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (26004,26005)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 0.279GB.  
  PID: 26004;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 26005;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	28329255	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (26092)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_PRO-seq_2_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_2_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_2_R1_processed.fastq
Analysis complete for H9_PRO-seq_2_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 0.279GB.  
  PID: 26092;	Command: fastqc;	Return code: 0;	Memory used: 0.177GB

> `FastQC report r1`	fastqc/H9_PRO-seq_2_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_dedup.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq` (26242)
<pre>
[INFO][0m 3047417 duplicated records removed
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 2.043GB.  
  PID: 26242;	Command: seqkit;	Return code: 0;	Memory used: 2.043GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq` (26585,26586)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 2.043GB.  
  PID: 26585;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 26586;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	45847892.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	28610139.0	PEPPRO	_RES_

> `Duplicate_reads`	3047417.0	PEPPRO	_RES_

> `Pct_reads_too_short`	49.4497	PEPPRO	_RES_

> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R2_cutadapt.txt` (26721)
<pre>
</pre>
Command completed. Elapsed time: 0:01:32. Running peak memory: 2.043GB.  
  PID: 26721;	Command: cutadapt;	Return code: 0;	Memory used: 0.274GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (26871,26872)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 2.043GB.  
  PID: 26871;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 26872;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	30855836	PEPPRO	_RES_

> `Trim_loss_rate`	73.33	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (26911)
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
<pre>
Started analysis of H9_PRO-seq_2_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_2_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_2_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_2_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 2.043GB.  
  PID: 26911;	Command: fastqc;	Return code: 0;	Memory used: 0.171GB

> `FastQC report r1`	fastqc/H9_PRO-seq_2_R2_trimmed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.histogram`  

> `fastq_pair -t 104142623 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq` (27149)
<pre>
Left paired: 29044502		Right paired: 29044502
Left single: 202372		Right single: 2678154
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:03:27. Running peak memory: 6.567GB.  
  PID: 27149;	Command: fastq_pair;	Return code: 0;	Memory used: 6.567GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq.paired.fq -o H9_PRO-seq_2 -d /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt` (27736)
<pre>
</pre>
Command completed. Elapsed time: 0:01:35. Running peak memory: 6.567GB.  
  PID: 27736;	Command: flash;	Return code: 0;	Memory used: 0.085GB


### Plot adapter insertion distribution (02-18 11:30:03) elapsed: 849.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt -u 8` (28233)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.567GB.  
  PID: 28233;	Command: Rscript;	Return code: 0;	Memory used: 0.238GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_2_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_2_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:30:08) elapsed: 5.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `Degradation_ratio`	1.0321	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq` (28260)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 6.567GB.  
  PID: 28260;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq.paired.fq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq.paired.fq`  

> `fastq_pair -t 104142623 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (28295)
<pre>
Left paired: 28105674		Right paired: 28105674
Left single: 223581		Right single: 2750162
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:34. Running peak memory: 6.567GB.  
  PID: 28295;	Command: fastq_pair;	Return code: 0;	Memory used: 5.373GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/repaired.flag`  

> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (28506)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 28506;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (28507)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 28507;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/repaired.flag` (28509)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 28509;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq.paired.fq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq`  

> `fastq_pair -t 104142623 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq` (28510)
<pre>
Left paired: 25286654		Right paired: 25286654
Left single: 163446		Right single: 5569182
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:05. Running peak memory: 6.567GB.  
  PID: 28510;	Command: fastq_pair;	Return code: 0;	Memory used: 5.733GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/dups_repaired.flag`  

> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq` (28863)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 28863;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq` (28865)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 28865;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/dups_repaired.flag` (28866)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 28866;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-18 11:35:13) elapsed: 305.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:35:13) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2` (28867)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 28867;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq` (28868)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_2 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
28105674 reads; of these:
  28105674 (100.00%) were unpaired; of these:
    25192086 (89.63%) aligned 0 times
    2913588 (10.37%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.37% overall alignment rate

> `Aligned_reads_human_rDNA`	5827176.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	18.89	PEPPRO	_RES_

### Map to human_rDNA (02-18 11:38:20) elapsed: 187.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2` (29086)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 29086;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R2.fq` (29087)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_2 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2913588 reads skipped
0 reads lost

### Map to genome (02-18 11:41:17) elapsed: 177.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_2 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/tmpeb_jq3k_ -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam` (29507,29508,29509)
<pre>
2252266 reads skipped
0 reads lost
25192086 reads; of these:
  25192086 (100.00%) were paired; of these:
    10443566 (41.46%) aligned concordantly 0 times
    12095912 (48.01%) aligned concordantly exactly 1 time
    2652608 (10.53%) aligned concordantly >1 times
    ----
    10443566 pairs aligned concordantly 0 times; of these:
      2323449 (22.25%) aligned discordantly 1 time
    ----
    8120117 pairs aligned 0 times concordantly or discordantly; of these:
      16240234 mates make up the pairs; of these:
        6664276 (41.04%) aligned 0 times
        3498753 (21.54%) aligned exactly 1 time
        6077205 (37.42%) aligned >1 times
86.77% overall alignment rate
[bam_sort_core] merging from 13 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 1:02:57. Running peak memory: 6.567GB.  
  PID: 29507;	Command: bowtie2;	Return code: 0;	Memory used: 3.676GB  
  PID: 29508;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 29509;	Command: samtools;	Return code: 0;	Memory used: 0.89GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (36213)
<pre>
</pre>
Command completed. Elapsed time: 0:01:30. Running peak memory: 6.567GB.  
  PID: 36213;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	43719896	PEPPRO	_RES_

> `QC_filtered_reads`	25554169	PEPPRO	_RES_

> `Aligned_reads`	18165727.0	PEPPRO	_RES_

> `Alignment_rate`	58.87	PEPPRO	_RES_

> `Total_efficiency`	15.7	PEPPRO	_RES_

> `Read_depth`	3.11	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_2 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/tmpeb_jq3k_ -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam` (37759,37760,37761)
<pre>
23034388 reads; of these:
  23034388 (100.00%) were paired; of these:
    9342849 (40.56%) aligned concordantly 0 times
    11259901 (48.88%) aligned concordantly exactly 1 time
    2431638 (10.56%) aligned concordantly >1 times
    ----
    9342849 pairs aligned concordantly 0 times; of these:
      2170368 (23.23%) aligned discordantly 1 time
    ----
    7172481 pairs aligned 0 times concordantly or discordantly; of these:
      14344962 mates make up the pairs; of these:
        6003927 (41.85%) aligned 0 times
        3261117 (22.73%) aligned exactly 1 time
        5079918 (35.41%) aligned >1 times
86.97% overall alignment rate
[bam_sort_core] merging from 12 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:56:15. Running peak memory: 6.567GB.  
  PID: 37759;	Command: bowtie2;	Return code: 0;	Memory used: 3.662GB  
  PID: 37760;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 37761;	Command: samtools;	Return code: 0;	Memory used: 0.873GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam` (45345)
<pre>
</pre>
Command completed. Elapsed time: 0:01:24. Running peak memory: 6.567GB.  
  PID: 45345;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (02-18 13:56:03) elapsed: 8086.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq` (45664)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 6.567GB.  
  PID: 45664;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq` (45715)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.567GB.  
  PID: 45715;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam` (45753)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 6.567GB.  
  PID: 45753;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	840435	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (45812)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.567GB.  
  PID: 45812;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_noMT.bam` (45850,45852,45853,45854)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 6.567GB.  
  PID: 45853;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 45850;	Command: samtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 45852;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 45854;	Command: xargs;	Return code: 0;	Memory used: 0.066GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (45903)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 45903;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (45904)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.567GB.  
  PID: 45904;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


### Split BAM file (02-18 13:59:14) elapsed: 191.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam` (45931,45932)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:33. Running peak memory: 6.567GB.  
  PID: 45931;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 45932;	Command: samtools;	Return code: 0;	Memory used: 4.687GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE2.bam` (46405,46407)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:11. Running peak memory: 6.567GB.  
  PID: 46405;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 46407;	Command: samtools;	Return code: 0;	Memory used: 3.758GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam` (46672)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 6.567GB.  
  PID: 46672;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam` (46900,46901)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:27. Running peak memory: 6.567GB.  
  PID: 46901;	Command: samtools;	Return code: 0;	Memory used: 4.083GB  
  PID: 46900;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE2.bam` (47100,47101)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:04. Running peak memory: 6.567GB.  
  PID: 47100;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 47101;	Command: samtools;	Return code: 0;	Memory used: 3.566GB


### Calculate library complexity (02-18 14:09:43) elapsed: 629.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_out.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam` (47254)
<pre>
BAM_INPUT
TOTAL READS     = 18471612
COUNTS_SUM      = 18471612
DISTINCT READS  = 1.4838e+07
DISTINCT COUNTS = 324
MAX COUNT       = 27760
COUNTS OF 1     = 1.31407e+07
OBSERVED COUNTS (27761)
1	13140744
2	1146563
3	271114
4	106759
5	54977
6	32098
7	20377
8	13755
9	9695
10	7243
11	5405
12	4178
13	3366
14	2630
15	2184
16	1670
17	1556
18	1258
19	1062
20	922
21	848
22	719
23	630
24	559
25	515
26	451
27	428
28	367
29	288
30	299
31	270
32	267
33	217
34	243
35	187
36	184
37	185
38	182
39	115
40	118
41	128
42	138
43	109
44	108
45	116
46	94
47	93
48	89
49	86
50	74
51	68
52	83
53	92
54	70
55	66
56	62
57	55
58	57
59	51
60	52
61	61
62	42
63	45
64	39
65	35
66	29
67	50
68	36
69	25
70	39
71	35
72	28
73	31
74	32
75	27
76	27
77	24
78	27
79	28
80	21
81	23
82	16
83	21
84	29
85	23
86	16
87	19
88	19
89	17
90	12
91	10
92	17
93	11
94	17
95	10
96	17
97	22
98	19
99	17
100	9
101	14
102	19
103	10
104	8
105	9
106	9
107	16
108	8
109	8
110	5
111	8
112	5
113	13
114	12
115	8
116	12
117	7
118	4
119	9
120	7
121	6
122	5
123	11
124	8
125	2
126	5
127	5
128	7
129	6
130	4
131	5
132	6
133	5
134	4
135	2
136	3
137	4
138	11
139	5
140	9
141	3
142	2
143	5
144	5
145	4
146	2
147	6
148	2
149	5
150	5
151	2
152	4
153	3
154	3
155	4
156	5
157	4
158	5
159	1
160	3
161	4
162	6
163	5
164	2
165	5
166	6
167	1
168	1
169	2
170	2
171	3
172	3
173	1
174	3
175	2
176	1
177	5
178	2
179	3
180	2
181	4
183	1
184	2
186	3
187	1
188	1
190	3
192	4
194	4
196	3
197	2
198	3
199	2
200	3
201	2
202	1
203	1
204	4
206	1
207	4
208	4
210	3
212	3
213	3
215	2
219	1
220	3
221	2
222	2
223	1
224	2
225	3
226	3
228	1
230	1
232	1
233	2
234	1
235	2
236	2
238	2
240	1
242	2
244	2
246	1
249	1
253	2
255	1
256	1
257	1
258	2
259	2
261	1
263	2
265	1
268	1
269	2
270	1
273	1
277	1
279	1
282	1
283	1
284	1
290	1
292	1
299	1
302	2
303	1
306	1
307	1
310	1
311	1
312	1
314	1
315	1
317	1
318	1
319	2
320	1
321	1
322	1
325	1
328	2
331	1
335	2
338	1
341	1
345	1
350	1
355	1
360	2
370	1
373	1
374	1
390	2
391	1
394	1
397	1
412	1
415	1
426	1
429	1
442	1
444	1
452	1
488	1
494	2
502	1
505	1
509	1
515	1
549	2
551	1
584	1
694	1
701	1
703	1
708	1
710	1
712	1
714	1
809	1
881	1
918	1
962	1
1001	1
1158	1
1264	1
1324	1
1628	1
1848	1
1890	1
1899	1
2227	1
2457	1
2584	1
3405	1
4384	1
9005	1
10491	1
13480	1
23174	1
27760	1

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
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 6.567GB.  
  PID: 47254;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam` (47644)
<pre>
BAM_INPUT
TOTAL READS     = 18471612
DISTINCT READS  = 1.4838e+07
DISTINCT COUNTS = 324
MAX COUNT       = 27760
COUNTS OF 1     = 1.31407e+07
MAX TERMS       = 100
OBSERVED COUNTS (27761)
1	13140744
2	1146563
3	271114
4	106759
5	54977
6	32098
7	20377
8	13755
9	9695
10	7243
11	5405
12	4178
13	3366
14	2630
15	2184
16	1670
17	1556
18	1258
19	1062
20	922
21	848
22	719
23	630
24	559
25	515
26	451
27	428
28	367
29	288
30	299
31	270
32	267
33	217
34	243
35	187
36	184
37	185
38	182
39	115
40	118
41	128
42	138
43	109
44	108
45	116
46	94
47	93
48	89
49	86
50	74
51	68
52	83
53	92
54	70
55	66
56	62
57	55
58	57
59	51
60	52
61	61
62	42
63	45
64	39
65	35
66	29
67	50
68	36
69	25
70	39
71	35
72	28
73	31
74	32
75	27
76	27
77	24
78	27
79	28
80	21
81	23
82	16
83	21
84	29
85	23
86	16
87	19
88	19
89	17
90	12
91	10
92	17
93	11
94	17
95	10
96	17
97	22
98	19
99	17
100	9
101	14
102	19
103	10
104	8
105	9
106	9
107	16
108	8
109	8
110	5
111	8
112	5
113	13
114	12
115	8
116	12
117	7
118	4
119	9
120	7
121	6
122	5
123	11
124	8
125	2
126	5
127	5
128	7
129	6
130	4
131	5
132	6
133	5
134	4
135	2
136	3
137	4
138	11
139	5
140	9
141	3
142	2
143	5
144	5
145	4
146	2
147	6
148	2
149	5
150	5
151	2
152	4
153	3
154	3
155	4
156	5
157	4
158	5
159	1
160	3
161	4
162	6
163	5
164	2
165	5
166	6
167	1
168	1
169	2
170	2
171	3
172	3
173	1
174	3
175	2
176	1
177	5
178	2
179	3
180	2
181	4
183	1
184	2
186	3
187	1
188	1
190	3
192	4
194	4
196	3
197	2
198	3
199	2
200	3
201	2
202	1
203	1
204	4
206	1
207	4
208	4
210	3
212	3
213	3
215	2
219	1
220	3
221	2
222	2
223	1
224	2
225	3
226	3
228	1
230	1
232	1
233	2
234	1
235	2
236	2
238	2
240	1
242	2
244	2
246	1
249	1
253	2
255	1
256	1
257	1
258	2
259	2
261	1
263	2
265	1
268	1
269	2
270	1
273	1
277	1
279	1
282	1
283	1
284	1
290	1
292	1
299	1
302	2
303	1
306	1
307	1
310	1
311	1
312	1
314	1
315	1
317	1
318	1
319	2
320	1
321	1
322	1
325	1
328	2
331	1
335	2
338	1
341	1
345	1
350	1
355	1
360	2
370	1
373	1
374	1
390	2
391	1
394	1
397	1
412	1
415	1
426	1
429	1
442	1
444	1
452	1
488	1
494	2
502	1
505	1
509	1
515	1
549	2
551	1
584	1
694	1
701	1
703	1
708	1
710	1
712	1
714	1
809	1
881	1
918	1
962	1
1001	1
1158	1
1264	1
1324	1
1628	1
1848	1
1890	1
1899	1
2227	1
2457	1
2584	1
3405	1
4384	1
9005	1
10491	1
13480	1
23174	1
27760	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
............................_.............................._............_.........................._....
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:52. Running peak memory: 6.567GB.  
  PID: 47644;	Command: preseq;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_counts.txt` (47761)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.567GB.  
  PID: 47761;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_plot` (47847)
<pre>
Processing H9_PRO-seq_2
INFO: Found real counts for H9_PRO-seq_2 - Total (M): 19.4822 Unique (M): 18.471612

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 6.567GB.  
  PID: 47847;	Command: Rscript;	Return code: 0;	Memory used: 0.371GB

> `Library complexity`	QC_hg38/H9_PRO-seq_2_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_2_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.854	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 14:13:50) elapsed: 247.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam` (47887)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.567GB.  
  PID: 47887;	Command: samtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv` (47903)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/tmp_H9_PRO-seq_2_PE1_9x4j0gv8'
Processing with 8 cores...
Discarding 100 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270755v1']
Keeping 95 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 6.567GB.  
  PID: 47903;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.097GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	9741100.0	PEPPRO	_RES_

> `PBC2`	9741100.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_unmap.bam` (47948)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.567GB.  
  PID: 47948;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam`

> `Unmapped_reads`	6664276	PEPPRO	_RES_

### Split BAM by strand (02-18 14:14:42) elapsed: 52.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam` (47991)
<pre>
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 6.567GB.  
  PID: 47991;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam` (48301)
<pre>
</pre>
Command completed. Elapsed time: 0:01:06. Running peak memory: 6.567GB.  
  PID: 48301;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 14:16:56) elapsed: 134.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (48416)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 48416;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt` (48417)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.567GB.  
  PID: 48417;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.503GB


> `TSS_Plus_Score`	33.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt` (48452)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.567GB.  
  PID: 48452;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.498GB


> `TSS_Minus_Score`	12.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt` (48480)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.567GB.  
  PID: 48480;	Command: Rscript;	Return code: 0;	Memory used: 0.236GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_2_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_2_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt` (48499,48500,48501,48502)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 48499;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 48501;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 48500;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 48502;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt` (48504)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 48504;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (02-18 14:17:17) elapsed: 21.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed` (48506,48507)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.567GB.  
  PID: 48506;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 48507;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed` (48510,48511)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 48510;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 48511;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSS_density.bed` (48513,48514,48515,48516)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.567GB.  
  PID: 48514;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 48516;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 48513;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 48515;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_gene_body_density.bed` (48557,48558,48559)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 6.567GB.  
  PID: 48559;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 48557;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB  
  PID: 48558;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed` (48600,48602,48603)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 48600;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 48603;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 48602;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	41.45	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed` (48608)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.567GB.  
  PID: 48608;	Command: Rscript;	Return code: 0;	Memory used: 0.237GB

> `Pause index`	QC_hg38/H9_PRO-seq_2_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_2_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed` (48628)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 48628;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 14:18:25) elapsed: 68.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam`
18165727.0 7048811

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam`
18165727.0 6572286

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_gene_sort.bed` (48702,48703)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 48702;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 48703;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_gene_coverage.bed` (48706)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 6.567GB.  
  PID: 48706;	Command: bedtools;	Return code: 0;	Memory used: 0.066GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed.gz` (48769)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 48769;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed` (48770)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 48770;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 14:19:34) elapsed: 69.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed` (48779)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 48779;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR"` (48781)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 48781;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed` (48782,48783,48784,48785)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 48782;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 48783;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 48785;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 48784;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed` (48787)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 48787;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_minus_coverage.bed` (48800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.567GB.  
  PID: 48800;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR"` (48844)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 48844;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed` (48857,48858,48861,48862)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 48857;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 48858;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 48862;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 48861;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed` (48955)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 48955;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_minus_coverage.bed` (49069)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.567GB.  
  PID: 49069;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed` (49090,49091,49092,49093)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 49090;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 49091;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 49093;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 49092;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed` (49098)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 49098;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_minus_coverage.bed` (49119)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.567GB.  
  PID: 49119;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed` (49131,49132,49133,49134)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.567GB.  
  PID: 49131;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 49132;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 49134;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 49133;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed` (49139)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.567GB.  
  PID: 49139;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_minus_coverage.bed` (49158)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 49158;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed` (49171,49172,49173,49174)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 49171;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 49173;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 49172;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 49174;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed` (49178)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.567GB.  
  PID: 49178;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_minus_coverage.bed` (49239)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.567GB.  
  PID: 49239;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed` (49288,49289,49290,49291)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 49288;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 49289;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 49291;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 49290;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed` (49293)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 49293;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_minus_coverage.bed` (49309)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 49309;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region"` (49330)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 49330;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed` (49331,49332,49333,49334)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 49331;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 49333;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 49332;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 49334;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (49337)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 49337;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed` (49358)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.567GB.  
  PID: 49358;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


### Plot FRiF/PRiF (02-18 14:22:49) elapsed: 195.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_2 -z 3099922541 -n 9944138 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (49390)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.567GB.  
  PID: 49390;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `FRiF`	QC_hg38/H9_PRO-seq_2_frif.pdf	FRiF	QC_hg38/H9_PRO-seq_2_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_2 -z 3099922541 -n 9944138 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (49450)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.567GB.  
  PID: 49450;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `PRiF`	QC_hg38/H9_PRO-seq_2_prif.pdf	PRiF	QC_hg38/H9_PRO-seq_2_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 14:23:47) elapsed: 58.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_exons_sort.bed` (49498,49499)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.567GB.  
  PID: 49499;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB  
  PID: 49498;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_introns_sort.bed` (49505,49506,49507)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.567GB.  
  PID: 49505;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 49507;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 49506;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_coverage.bed` (49513)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 6.567GB.  
  PID: 49513;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_coverage.bed` (49557)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.567GB.  
  PID: 49557;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.165727)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_rpkm.bed` (49598,49599,49600)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 49598;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 49600;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 49599;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.165727)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_rpkm.bed` (49603,49604,49605)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 49603;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 49605;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 49604;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed` (49608,49609,49610)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 49608;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 49610;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 49609;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed --annotate` (49617)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 6.567GB.  
  PID: 49617;	Command: Rscript;	Return code: 0;	Memory used: 0.366GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_2_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_2_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed` (49644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 49644;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-18 14:24:55) elapsed: 69.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam` (49652)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.567GB.  
  PID: 49652;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (49863)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_2_plus_cuttrace_myp4uy5k'
Processing with 2 cores...
stdin is empty of data
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:55. Running peak memory: 6.567GB.  
  PID: 49863;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.603GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam` (52622)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.567GB.  
  PID: 52622;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_exact_body_0-mer.bw -p 5 --variable-step --tail-edge` (52824)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_2_minus_cuttrace_e7oy3pli'
Processing with 5 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr14_KI270723v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270755v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.567GB.  
  PID: 52824;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.33GB

Starting cleanup: 68 files; 4 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:26:25
*  Total elapsed time (all runs):  6:56:14
*         Peak memory (this run):  6.5669 GB
*        Pipeline completed time: 2020-02-18 14:35:23
