### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name H9_treated_PRO-seq_1 --genome hg38 --input /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz --protocol PRO --umi-len 8 --max-len -1 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj40-16c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/
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
*              `input`:  `['/project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_treated_PRO-seq_1`
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

Local input file: /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz
Local input file: /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz

> `File_mb`	2710.39	PEPPRO	_RES_

> `Read_type`	paired	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-18 11:08:59) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz` (357287)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 357287;	Command: ln;	Return code: 0;	Memory used: 0.002GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz` (357290)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 357290;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1.fastq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1.fastq` (357293)
<pre>
</pre>
Command completed. Elapsed time: 0:02:16. Running peak memory: 0.003GB.  
  PID: 357293;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2.fastq` (231032)
<pre>
</pre>
Command completed. Elapsed time: 0:01:40. Running peak memory: 0.003GB.  
  PID: 231032;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	111300168	PEPPRO	_RES_

> `Fastq_reads`	111300168	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz']

### FASTQ processing:  (02-18 11:15:10) elapsed: 372.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt` (221538)
<pre>
</pre>
Command completed. Elapsed time: 0:01:40. Running peak memory: 0.288GB.  
  PID: 221538;	Command: cutadapt;	Return code: 0;	Memory used: 0.288GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (221693,221694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 0.288GB.  
  PID: 221693;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 221694;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	29119715	PEPPRO	_RES_

> `Trim_loss_rate`	73.84	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (221768)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_treated_PRO-seq_1_R1_processed.fastq
Approx 5% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 10% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 15% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 20% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 25% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 30% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 35% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 40% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 45% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 50% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 55% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 60% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 65% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 70% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 75% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 80% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 85% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 90% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Approx 95% complete for H9_treated_PRO-seq_1_R1_processed.fastq
Analysis complete for H9_treated_PRO-seq_1_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 0.288GB.  
  PID: 221768;	Command: fastqc;	Return code: 0;	Memory used: 0.175GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_dedup.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq` (221916)
<pre>
[INFO][0m 2891628 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 2.591GB.  
  PID: 221916;	Command: seqkit;	Return code: 0;	Memory used: 2.591GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq` (222009,222010)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 2.591GB.  
  PID: 222009;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 222010;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	43631018.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_too_short`	25201427.0	PEPPRO	_RES_

> `Duplicate_reads`	2891628.0	PEPPRO	_RES_

> `Pct_reads_too_short`	45.2855	PEPPRO	_RES_

> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2.fastq -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R2_cutadapt.txt` (222387)
<pre>
</pre>
Command completed. Elapsed time: 0:01:30. Running peak memory: 2.591GB.  
  PID: 222387;	Command: cutadapt;	Return code: 0;	Memory used: 0.273GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (222530,222531)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 2.591GB.  
  PID: 222530;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 222531;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	31193330	PEPPRO	_RES_

> `Trim_loss_rate`	71.97	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (222576)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 5% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 10% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 15% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 20% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 25% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 30% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 35% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 40% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 45% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 50% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 55% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 60% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 65% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 70% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 75% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 80% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 85% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 90% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Approx 95% complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
Analysis complete for H9_treated_PRO-seq_1_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 2.591GB.  
  PID: 222576;	Command: fastqc;	Return code: 0;	Memory used: 0.171GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_1_R2_trimmed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.histogram`  

> `fastq_pair -t 100170151 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq` (222683)
<pre>
Left paired: 30307647		Right paired: 30307647
Left single: 141010		Right single: 2158148
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:57. Running peak memory: 6.727GB.  
  PID: 222683;	Command: fastq_pair;	Return code: 0;	Memory used: 6.727GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq.paired.fq -o H9_treated_PRO-seq_1 -d /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt` (223219)
<pre>
</pre>
Command completed. Elapsed time: 0:01:35. Running peak memory: 6.727GB.  
  PID: 223219;	Command: flash;	Return code: 0;	Memory used: 0.125GB


### Plot adapter insertion distribution (02-18 11:29:02) elapsed: 832.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt -u 8` (223464)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.727GB.  
  PID: 223464;	Command: Rscript;	Return code: 0;	Memory used: 0.127GB

> `Adapter insertion distribution`	cutadapt/H9_treated_PRO-seq_1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_treated_PRO-seq_1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-18 11:29:10) elapsed: 8.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `Degradation_ratio`	1.158	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq` (223511)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.727GB.  
  PID: 223511;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq.paired.fq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq.paired.fq`  

> `fastq_pair -t 100170151 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (223537)
<pre>
Left paired: 28948010		Right paired: 28948010
Left single: 171705		Right single: 2245320
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:25. Running peak memory: 6.727GB.  
  PID: 223537;	Command: fastq_pair;	Return code: 0;	Memory used: 6.024GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/repaired.flag`  

> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (224083)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 224083;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (224084)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 224084;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/repaired.flag` (224086)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 224086;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq.paired.fq`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq`  

> `fastq_pair -t 100170151 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq` (224087)
<pre>
Left paired: 26293117		Right paired: 26293117
Left single: 111289		Right single: 4900213
Writing the paired reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:02. Running peak memory: 6.727GB.  
  PID: 224087;	Command: fastq_pair;	Return code: 0;	Memory used: 6.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/dups_repaired.flag`  

> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq` (224237)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 224237;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq` (224241)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 224241;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/dups_repaired.flag` (224243)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 224243;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-18 11:34:09) elapsed: 299.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-18 11:34:09) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2` (224246)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 224246;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq` (224247)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
28948010 reads; of these:
  28948010 (100.00%) were unpaired; of these:
    25898407 (89.47%) aligned 0 times
    3049603 (10.53%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.53% overall alignment rate

> `Aligned_reads_human_rDNA`	6099206.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	19.55	PEPPRO	_RES_

### Map to human_rDNA (02-18 11:37:30) elapsed: 201.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2` (224672)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 224672;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R2.fq` (224673)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
3049603 reads skipped
0 reads lost

### Map to genome (02-18 11:40:35) elapsed: 185.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_1 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/tmplu04a8tr -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam` (225114,225115,225117)
<pre>
2503948 reads skipped
0 reads lost
25898407 reads; of these:
  25898407 (100.00%) were paired; of these:
    12391166 (47.85%) aligned concordantly 0 times
    11262767 (43.49%) aligned concordantly exactly 1 time
    2244474 (8.67%) aligned concordantly >1 times
    ----
    12391166 pairs aligned concordantly 0 times; of these:
      2903829 (23.43%) aligned discordantly 1 time
    ----
    9487337 pairs aligned 0 times concordantly or discordantly; of these:
      18974674 mates make up the pairs; of these:
        7715519 (40.66%) aligned 0 times
        3997151 (21.07%) aligned exactly 1 time
        7262004 (38.27%) aligned >1 times
85.10% overall alignment rate
[bam_sort_core] merging from 14 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 1:01:45. Running peak memory: 6.727GB.  
  PID: 225114;	Command: bowtie2;	Return code: 0;	Memory used: 3.693GB  
  PID: 225115;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 225117;	Command: samtools;	Return code: 0;	Memory used: 0.875GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_fail_qc.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (231266)
<pre>
</pre>
Command completed. Elapsed time: 0:01:39. Running peak memory: 6.727GB.  
  PID: 231266;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	44081295	PEPPRO	_RES_

> `QC_filtered_reads`	26072183	PEPPRO	_RES_

> `Aligned_reads`	18009111.5	PEPPRO	_RES_

> `Alignment_rate`	57.73	PEPPRO	_RES_

> `Total_efficiency`	16.18	PEPPRO	_RES_

> `Read_depth`	3.41	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_1 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/tmplu04a8tr -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam` (232653,232654,232655)
<pre>
23789169 reads; of these:
  23789169 (100.00%) were paired; of these:
    11229616 (47.20%) aligned concordantly 0 times
    10491312 (44.10%) aligned concordantly exactly 1 time
    2068241 (8.69%) aligned concordantly >1 times
    ----
    11229616 pairs aligned concordantly 0 times; of these:
      2715002 (24.18%) aligned discordantly 1 time
    ----
    8514614 pairs aligned 0 times concordantly or discordantly; of these:
      17029228 mates make up the pairs; of these:
        7039185 (41.34%) aligned 0 times
        3733545 (21.92%) aligned exactly 1 time
        6256498 (36.74%) aligned >1 times
85.21% overall alignment rate
[bam_sort_core] merging from 12 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:54:05. Running peak memory: 6.727GB.  
  PID: 232654;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 232653;	Command: bowtie2;	Return code: 0;	Memory used: 3.666GB  
  PID: 232655;	Command: samtools;	Return code: 0;	Memory used: 0.872GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam` (237866)
<pre>
</pre>
Command completed. Elapsed time: 0:01:36. Running peak memory: 6.727GB.  
  PID: 237866;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (02-18 13:51:21) elapsed: 7846.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq` (238192)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 6.727GB.  
  PID: 238192;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq` (238245)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.727GB.  
  PID: 238245;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam` (238283)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 6.727GB.  
  PID: 238283;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	655778	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (238328)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.727GB.  
  PID: 238328;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_noMT.bam` (238366,238367,238368,238370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 6.727GB.  
  PID: 238368;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 238366;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 238367;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 238370;	Command: xargs;	Return code: 0;	Memory used: 0.081GB


> `mv /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_noMT.bam /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (238422)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 238422;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (238423)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.727GB.  
  PID: 238423;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


### Split BAM file (02-18 13:54:42) elapsed: 200.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam` (238458,238459)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:35. Running peak memory: 6.727GB.  
  PID: 238458;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 238459;	Command: samtools;	Return code: 0;	Memory used: 4.884GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE2.bam` (238858,238859)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:06. Running peak memory: 6.727GB.  
  PID: 238858;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 238859;	Command: samtools;	Return code: 0;	Memory used: 3.788GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam` (239306)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 6.727GB.  
  PID: 239306;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam` (239342,239343)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:25. Running peak memory: 6.727GB.  
  PID: 239343;	Command: samtools;	Return code: 0;	Memory used: 4.3GB  
  PID: 239342;	Command: samtools;	Return code: 0;	Memory used: 0.004GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE2.bam` (239579,239580)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:00. Running peak memory: 6.727GB.  
  PID: 239579;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 239580;	Command: samtools;	Return code: 0;	Memory used: 3.582GB


### Calculate library complexity (02-18 14:05:02) elapsed: 620.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_out.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam` (239848)
<pre>
BAM_INPUT
TOTAL READS     = 18664481
COUNTS_SUM      = 18664481
DISTINCT READS  = 1.40562e+07
DISTINCT COUNTS = 291
MAX COUNT       = 20745
COUNTS OF 1     = 1.22495e+07
OBSERVED COUNTS (20746)
1	12249493
2	1104689
3	295074
4	133371
5	74603
6	47569
7	32572
8	23059
9	16807
10	12896
11	10229
12	7936
13	6444
14	5351
15	4479
16	3679
17	3062
18	2646
19	2181
20	1837
21	1659
22	1506
23	1264
24	1153
25	1027
26	944
27	768
28	727
29	604
30	557
31	516
32	468
33	453
34	389
35	387
36	322
37	310
38	267
39	272
40	254
41	240
42	242
43	218
44	185
45	176
46	145
47	151
48	133
49	124
50	125
51	122
52	98
53	98
54	88
55	88
56	88
57	77
58	88
59	74
60	79
61	68
62	66
63	49
64	54
65	57
66	40
67	48
68	34
69	48
70	43
71	38
72	35
73	27
74	25
75	38
76	35
77	24
78	32
79	30
80	36
81	24
82	22
83	31
84	16
85	32
86	20
87	14
88	13
89	12
90	8
91	15
92	6
93	13
94	17
95	19
96	15
97	10
98	9
99	14
100	17
101	13
102	8
103	7
104	14
105	9
106	8
107	7
108	11
109	12
110	7
111	10
112	12
113	5
114	12
115	7
116	11
117	4
118	7
119	8
120	3
121	10
122	11
123	8
124	3
125	8
126	6
127	6
128	2
129	6
130	3
131	3
132	5
133	10
134	2
135	4
136	5
137	4
138	2
139	2
140	2
141	3
142	8
143	3
144	2
145	3
146	4
147	5
148	3
149	3
150	5
152	1
154	2
155	2
156	3
157	4
158	2
159	5
160	3
161	2
162	3
163	5
164	5
166	3
167	6
168	4
169	3
170	1
171	3
172	7
173	2
174	1
175	3
176	1
177	3
178	2
179	1
180	1
181	1
182	2
183	1
184	2
185	3
186	1
187	1
188	1
189	2
190	2
191	2
192	2
194	1
195	1
196	1
197	1
198	2
199	1
200	2
201	2
203	1
204	1
205	1
206	1
208	3
210	1
212	1
213	1
214	1
215	6
217	1
219	1
220	1
221	2
222	2
226	3
229	1
232	3
235	1
238	1
239	1
240	2
241	2
248	2
249	1
254	2
256	1
259	1
261	1
262	1
264	1
265	1
268	2
271	1
272	1
273	1
275	1
285	1
286	1
295	1
296	2
299	1
300	1
301	1
303	1
308	1
313	2
314	1
317	1
334	1
339	1
347	1
348	1
350	1
354	1
361	1
376	1
382	1
391	1
402	1
408	1
411	1
426	1
445	1
451	1
474	1
502	1
506	1
513	1
598	1
614	1
644	1
653	1
710	1
763	1
766	1
796	1
1031	1
1161	1
1187	1
1340	1
1349	1
1417	1
1428	1
1705	1
1759	1
1901	1
2597	1
3806	1
7544	1
7762	1
9871	1
20534	1
20745	1

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
Command completed. Elapsed time: 0:01:42. Running peak memory: 6.727GB.  
  PID: 239848;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam` (240065)
<pre>
BAM_INPUT
TOTAL READS     = 18664481
DISTINCT READS  = 1.40562e+07
DISTINCT COUNTS = 291
MAX COUNT       = 20745
COUNTS OF 1     = 1.22495e+07
MAX TERMS       = 100
OBSERVED COUNTS (20746)
1	12249493
2	1104689
3	295074
4	133371
5	74603
6	47569
7	32572
8	23059
9	16807
10	12896
11	10229
12	7936
13	6444
14	5351
15	4479
16	3679
17	3062
18	2646
19	2181
20	1837
21	1659
22	1506
23	1264
24	1153
25	1027
26	944
27	768
28	727
29	604
30	557
31	516
32	468
33	453
34	389
35	387
36	322
37	310
38	267
39	272
40	254
41	240
42	242
43	218
44	185
45	176
46	145
47	151
48	133
49	124
50	125
51	122
52	98
53	98
54	88
55	88
56	88
57	77
58	88
59	74
60	79
61	68
62	66
63	49
64	54
65	57
66	40
67	48
68	34
69	48
70	43
71	38
72	35
73	27
74	25
75	38
76	35
77	24
78	32
79	30
80	36
81	24
82	22
83	31
84	16
85	32
86	20
87	14
88	13
89	12
90	8
91	15
92	6
93	13
94	17
95	19
96	15
97	10
98	9
99	14
100	17
101	13
102	8
103	7
104	14
105	9
106	8
107	7
108	11
109	12
110	7
111	10
112	12
113	5
114	12
115	7
116	11
117	4
118	7
119	8
120	3
121	10
122	11
123	8
124	3
125	8
126	6
127	6
128	2
129	6
130	3
131	3
132	5
133	10
134	2
135	4
136	5
137	4
138	2
139	2
140	2
141	3
142	8
143	3
144	2
145	3
146	4
147	5
148	3
149	3
150	5
152	1
154	2
155	2
156	3
157	4
158	2
159	5
160	3
161	2
162	3
163	5
164	5
166	3
167	6
168	4
169	3
170	1
171	3
172	7
173	2
174	1
175	3
176	1
177	3
178	2
179	1
180	1
181	1
182	2
183	1
184	2
185	3
186	1
187	1
188	1
189	2
190	2
191	2
192	2
194	1
195	1
196	1
197	1
198	2
199	1
200	2
201	2
203	1
204	1
205	1
206	1
208	3
210	1
212	1
213	1
214	1
215	6
217	1
219	1
220	1
221	2
222	2
226	3
229	1
232	3
235	1
238	1
239	1
240	2
241	2
248	2
249	1
254	2
256	1
259	1
261	1
262	1
264	1
265	1
268	2
271	1
272	1
273	1
275	1
285	1
286	1
295	1
296	2
299	1
300	1
301	1
303	1
308	1
313	2
314	1
317	1
334	1
339	1
347	1
348	1
350	1
354	1
361	1
376	1
382	1
391	1
402	1
408	1
411	1
426	1
445	1
451	1
474	1
502	1
506	1
513	1
598	1
614	1
644	1
653	1
710	1
763	1
766	1
796	1
1031	1
1161	1
1187	1
1340	1
1349	1
1417	1
1428	1
1705	1
1759	1
1901	1
2597	1
3806	1
7544	1
7762	1
9871	1
20534	1
20745	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
......................................._......_.__......_....................................._..._........
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:48. Running peak memory: 6.727GB.  
  PID: 240065;	Command: preseq;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam) > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_counts.txt` (240182)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.727GB.  
  PID: 240182;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_plot` (240225)
<pre>
Processing H9_treated_PRO-seq_1
INFO: Found real counts for H9_treated_PRO-seq_1 - Total (M): 19.73867 Unique (M): 18.664481

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.727GB.  
  PID: 240225;	Command: Rscript;	Return code: 0;	Memory used: 0.201GB

> `Library complexity`	QC_hg38/H9_treated_PRO-seq_1_preseq_plot.pdf	Library complexity	QC_hg38/H9_treated_PRO-seq_1_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8097	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-18 14:09:07) elapsed: 245.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam` (240258)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.727GB.  
  PID: 240258;	Command: samtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv` (240274)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/tmp_H9_treated_PRO-seq_1_PE1_pe0yqviy'
Processing with 8 cores...
Discarding 94 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000214v1']
Keeping 101 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.727GB.  
  PID: 240274;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.9GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	9869335.0	PEPPRO	_RES_

> `PBC2`	9869335.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_unmap.bam` (240339)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.727GB.  
  PID: 240339;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam`

> `Unmapped_reads`	7715519	PEPPRO	_RES_

### Split BAM by strand (02-18 14:10:04) elapsed: 57.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam` (240616)
<pre>
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 6.727GB.  
  PID: 240616;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam` (240708)
<pre>
</pre>
Command completed. Elapsed time: 0:01:07. Running peak memory: 6.727GB.  
  PID: 240708;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-18 14:12:19) elapsed: 135.0 _TIME_

Missing stat 'TSS_Minus_Score'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (240798)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 240798;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt` (240800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.727GB.  
  PID: 240800;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.614GB


> `TSS_Plus_Score`	53.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt` (240826)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.727GB.  
  PID: 240826;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.63GB


> `TSS_Minus_Score`	17.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt` (240851)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 6.727GB.  
  PID: 240851;	Command: Rscript;	Return code: 0;	Memory used: 0.346GB

> `TSS enrichment`	QC_hg38/H9_treated_PRO-seq_1_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_treated_PRO-seq_1_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt` (240871,240872,240873,240874)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 240871;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 240873;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 240872;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 240874;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt` (240876)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 240876;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-18 14:12:41) elapsed: 22.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed` (240878,240879)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.727GB.  
  PID: 240878;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 240879;	Command: bedtools;	Return code: 0;	Memory used: 0.1GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed` (240883,240884)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 240883;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 240884;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSS_density.bed` (240886,240887,240888,240889)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.727GB.  
  PID: 240887;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 240889;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 240886;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 240888;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_gene_body_density.bed` (240936,240937,240938)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 6.727GB.  
  PID: 240937;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 240936;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB  
  PID: 240938;	Command: sort;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSS_density.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed` (241016,241017,241018)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 241016;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 241018;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 241017;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	75.18	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed` (241024)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 6.727GB.  
  PID: 241024;	Command: Rscript;	Return code: 0;	Memory used: 0.367GB

> `Pause index`	QC_hg38/H9_treated_PRO-seq_1_pause_index.pdf	Pause index	QC_hg38/H9_treated_PRO-seq_1_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed` (241045)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 241045;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-18 14:13:46) elapsed: 65.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam`
18009111.5 7282779

> `Plus_FRiP`	0.4	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam`
18009111.5 6888359

> `Minus_FRiP`	0.38	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_gene_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_gene_sort.bed` (241115,241116)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 241115;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 241116;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_gene_coverage.bed` (241119)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.727GB.  
  PID: 241119;	Command: bedtools;	Return code: 0;	Memory used: 0.06GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed.gz` (241197)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 241197;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed` (241198)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 241198;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate fraction and proportion of reads in features (FRiF/PRiF) (02-18 14:14:54) elapsed: 68.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed` (241206)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 241206;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR"` (241208)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 241208;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed` (241209,241210,241211,241212)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 241209;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 241210;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 241212;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 241211;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed` (241215)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.727GB.  
  PID: 241215;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_minus_coverage.bed` (241432)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.727GB.  
  PID: 241432;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR"` (241448)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 241448;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed` (241449,241450,241451,241452)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 241449;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 241450;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 241452;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 241451;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed` (241455)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.727GB.  
  PID: 241455;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_minus_coverage.bed` (241467)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.727GB.  
  PID: 241467;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed` (241503,241504,241505,241506)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 241503;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 241504;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 241506;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB  
  PID: 241505;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed` (241509)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.727GB.  
  PID: 241509;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_minus_coverage.bed` (241523)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.727GB.  
  PID: 241523;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed` (241558,241559,241560,241561)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.727GB.  
  PID: 241558;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 241559;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 241561;	Command: bedtools;	Return code: 0;	Memory used: 0.174GB  
  PID: 241560;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed` (241566)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.727GB.  
  PID: 241566;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_minus_coverage.bed` (241601)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.727GB.  
  PID: 241601;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed` (241625,241626,241627,241628)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 241625;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 241627;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 241626;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 241628;	Command: bedtools;	Return code: 0;	Memory used: 0.085GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed` (241631)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.727GB.  
  PID: 241631;	Command: bedtools;	Return code: 0;	Memory used: 0.076GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_minus_coverage.bed` (241649)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.727GB.  
  PID: 241649;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed` (241670,241671,241672,241673)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 241670;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 241671;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 241673;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 241672;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed` (241675)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.727GB.  
  PID: 241675;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_minus_coverage.bed` (241690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.727GB.  
  PID: 241690;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target exists: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region"` (241711)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 241711;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed` (241712,241713,241714,241715)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 241712;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 241714;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 241713;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 241715;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (241718)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.727GB.  
  PID: 241718;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed` (241735)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.727GB.  
  PID: 241735;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


### Plot FRiF/PRiF (02-18 14:18:10) elapsed: 196.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_frif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_1 -z 3099922541 -n 10005525 -y frif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_frif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (241765)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.727GB.  
  PID: 241765;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `FRiF`	QC_hg38/H9_treated_PRO-seq_1_frif.pdf	FRiF	QC_hg38/H9_treated_PRO-seq_1_frif.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_prif.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_1 -z 3099922541 -n 10005525 -y prif --reads -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_prif.pdf --bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (241831)
<pre>
Cumulative prif plot completed!

</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.727GB.  
  PID: 241831;	Command: Rscript;	Return code: 0;	Memory used: 0.468GB

> `PRiF`	QC_hg38/H9_treated_PRO-seq_1_prif.pdf	PRiF	QC_hg38/H9_treated_PRO-seq_1_prif.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-18 14:19:07) elapsed: 57.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_exons_sort.bed` (241910,241911)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.727GB.  
  PID: 241911;	Command: bedtools;	Return code: 0;	Memory used: 0.092GB  
  PID: 241910;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_introns_sort.bed` (241917,241918,241919)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.727GB.  
  PID: 241917;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 241919;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 241918;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_coverage.bed` (241926)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.727GB.  
  PID: 241926;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_coverage.bed` (241963)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.727GB.  
  PID: 241963;	Command: bedtools;	Return code: 0;	Memory used: 0.075GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.0091115)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_rpkm.bed` (242221,242222,242223)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 242221;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 242223;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 242222;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.0091115)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_rpkm.bed` (242225,242226,242228)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 242225;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 242228;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 242226;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_rpkm.bed /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed` (242230,242231,242232)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 242230;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 242232;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 242231;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.19	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed --annotate` (242238)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.727GB.  
  PID: 242238;	Command: Rscript;	Return code: 0;	Memory used: 0.236GB

> `mRNA contamination`	QC_hg38/H9_treated_PRO-seq_1_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_treated_PRO-seq_1_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed` (242257)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 242257;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-18 14:20:18) elapsed: 71.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam` (242265)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.727GB.  
  PID: 242265;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (242272)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_1_plus_cuttrace_ykogbwot'
Processing with 2 cores...
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270752v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:46. Running peak memory: 6.727GB.  
  PID: 242272;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.607GB

Target to produce: `/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam` (245023)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.727GB.  
  PID: 245023;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_exact_body_0-mer.bw -p 5 --variable-step --tail-edge` (245030)
<pre>
Registering input file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_1_minus_cuttrace_8p7wjpdx'
Processing with 5 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 106 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr17_KI270730v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1']
Keeping 89 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 89 files into output file: '/project/shefflab/processed/peppro/paper/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_exact_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.727GB.  
  PID: 245030;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 0.32GB

Starting cleanup: 68 files; 4 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:21:40
*  Total elapsed time (all runs):  6:43:56
*         Peak memory (this run):  6.7273 GB
*        Pipeline completed time: 2020-02-18 14:30:38
