### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name H9_PRO-seq_2 --genome hg38 --input /project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz --protocol PRO --umi-len 8 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj37-15c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/
*  Pipeline started at:   (02-27 09:25:33) elapsed: 1.0 _TIME_

### Version log:

*       Python version:  3.6.5
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.1
*        Pipeline hash:  2fe0657f50e41000560af043f4914b3a240296f2
*      Pipeline branch:  * dev
*        Pipeline date:  2020-02-27 09:20:39 -0500

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
*      `output_parent`:  `/project/shefflab/processed/peppro/paper/dev4/results_pipeline`
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

### Merge/link and fastq conversion:  (02-27 09:25:33) elapsed: 0.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE1.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz` (308236)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 308236;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_DMSO_rep2_PE2.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz` (308237)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 308237;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq` (308238)
<pre>
</pre>
Command completed. Elapsed time: 0:01:05. Running peak memory: 0.003GB.  
  PID: 308238;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2.fastq` (122142)
<pre>
</pre>
Command completed. Elapsed time: 0:01:03. Running peak memory: 0.003GB.  
  PID: 122142;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	115714026	PEPPRO	_RES_

> `Fastq_reads`	115714026	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/H9_PRO-seq_2_R2.fastq.gz']

### FASTQ processing:  (02-27 09:29:20) elapsed: 226.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt` (422684)
<pre>
</pre>
Command completed. Elapsed time: 0:01:23. Running peak memory: 0.271GB.  
  PID: 422684;	Command: cutadapt;	Return code: 0;	Memory used: 0.271GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (423011,423012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 0.271GB.  
  PID: 423011;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 423012;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	28329255	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (423070)
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
Command completed. Elapsed time: 0:01:17. Running peak memory: 0.271GB.  
  PID: 423070;	Command: fastqc;	Return code: 0;	Memory used: 0.176GB

> `FastQC report r1`	fastqc/H9_PRO-seq_2_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_dedup.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq` (423162)
<pre>
[INFO][0m 3047417 duplicated records removed
</pre>
Command completed. Elapsed time: 0:01:11. Running peak memory: 2.039GB.  
  PID: 423162;	Command: seqkit;	Return code: 0;	Memory used: 2.039GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq` (423250,423251)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 2.039GB.  
  PID: 423250;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 423251;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	45847892.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	28610139.0	PEPPRO	_RES_

> `Duplicate_reads`	3047417.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.4497	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/processed_R1.flag` (423573)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.039GB.  
  PID: 423573;	Command: touch;	Return code: 0;	Memory used: 0.002GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_R2_cutadapt.txt` (423575)
<pre>
</pre>
Command completed. Elapsed time: 0:01:24. Running peak memory: 2.039GB.  
  PID: 423575;	Command: cutadapt;	Return code: 0;	Memory used: 0.254GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (423681,423682)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 2.039GB.  
  PID: 423681;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 423682;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	56658510	PEPPRO	_RES_

> `Trim_loss_rate`	51.04	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (423724)
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
Command completed. Elapsed time: 0:01:22. Running peak memory: 2.039GB.  
  PID: 423724;	Command: fastqc;	Return code: 0;	Memory used: 0.174GB

> `FastQC report r1`	fastqc/H9_PRO-seq_2_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (423827)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
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
Command completed. Elapsed time: 0:01:12. Running peak memory: 2.039GB.  
  PID: 423827;	Command: fastqc;	Return code: 0;	Memory used: 0.171GB

> `FastQC report r2`	fastqc/H9_PRO-seq_2_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.histogram`  

> `fastq_pair -t 104142623 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq` (423915)
<pre>
Left paired: 29044502		Right paired: 29044502
Left single: 202372		Right single: 2678154
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:04:06. Running peak memory: 6.567GB.  
  PID: 423915;	Command: fastq_pair;	Return code: 0;	Memory used: 6.567GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_noadap.fastq.paired.fq -o H9_PRO-seq_2 -d /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt` (424376)
<pre>
</pre>
Command completed. Elapsed time: 0:01:33. Running peak memory: 6.567GB.  
  PID: 424376;	Command: flash;	Return code: 0;	Memory used: 0.088GB


### Plot adapter insertion distribution (02-27 09:45:22) elapsed: 962.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt -u 8` (424816)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.567GB.  
  PID: 424816;	Command: Rscript;	Return code: 0;	Memory used: 0.124GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_2_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_2_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:45:29) elapsed: 8.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/cutadapt/H9_PRO-seq_2.hist`

> `Degradation_ratio`	1.0321	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq` (424845)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 6.567GB.  
  PID: 424845;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/processed_R2.flag` (424867)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 424867;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/repaired.flag`  

> `fastq_pair -t 104142623 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (424869)
<pre>
Left paired: 28105674		Right paired: 28105674
Left single: 223581		Right single: 2750162
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:32. Running peak memory: 6.567GB.  
  PID: 424869;	Command: fastq_pair;	Return code: 0;	Memory used: 5.726GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq` (425024)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 425024;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq` (425025)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 425025;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/repaired.flag` (425027)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 425027;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/dups_repaired.flag`  

> `fastq_pair -t 104142623 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq` (425028)
<pre>
Left paired: 25286654		Right paired: 25286654
Left single: 163446		Right single: 5569182
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:13. Running peak memory: 6.567GB.  
  PID: 425028;	Command: fastq_pair;	Return code: 0;	Memory used: 5.798GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq` (425376)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 425376;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq` (425377)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 425377;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/dups_repaired.flag` (425379)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 425379;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-27 09:50:40) elapsed: 310.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:50:40) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2` (425380)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 425380;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq` (425381)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_2 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
28105674 reads; of these:
  28105674 (100.00%) were unpaired; of these:
    25192086 (89.63%) aligned 0 times
    2913588 (10.37%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.37% overall alignment rate

> `Aligned_reads_human_rDNA`	5827176.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.28	PEPPRO	_RES_

### Map to human_rDNA (02-27 09:55:27) elapsed: 288.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2` (425879)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 425879;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R2.fq` (425880)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_2 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/fastq/H9_PRO-seq_2_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2913588 reads skipped
0 reads lost

### Map to genome (02-27 10:00:01) elapsed: 274.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_2 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/tmpz4pucrzs -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam` (426323,426336,426337)
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
Command completed. Elapsed time: 1:03:58. Running peak memory: 6.567GB.  
  PID: 426323;	Command: bowtie2;	Return code: 0;	Memory used: 3.673GB  
  PID: 426336;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 426337;	Command: samtools;	Return code: 0;	Memory used: 0.874GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (432870)
<pre>
</pre>
Command completed. Elapsed time: 0:01:41. Running peak memory: 6.567GB.  
  PID: 432870;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	43719896	PEPPRO	_RES_

> `QC_filtered_reads`	25554169	PEPPRO	_RES_

> `Aligned_reads`	18165727.0	PEPPRO	_RES_

> `Alignment_rate`	32.06	PEPPRO	_RES_

> `Total_efficiency`	15.7	PEPPRO	_RES_

> `Read_depth`	3.11	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_PRO-seq_2 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/tmpz4pucrzs -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam` (434438,434445,434446)
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
Command completed. Elapsed time: 0:57:56. Running peak memory: 6.567GB.  
  PID: 434438;	Command: bowtie2;	Return code: 0;	Memory used: 3.662GB  
  PID: 434445;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 434446;	Command: samtools;	Return code: 0;	Memory used: 0.884GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam` (440535)
<pre>
</pre>
Command completed. Elapsed time: 0:01:32. Running peak memory: 6.567GB.  
  PID: 440535;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (02-27 12:18:23) elapsed: 8301.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R2.fq` (440648)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.567GB.  
  PID: 440648;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/prealignments/H9_PRO-seq_2_human_rDNA_unmap_R1.fq` (440690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 6.567GB.  
  PID: 440690;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam` (440735)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 6.567GB.  
  PID: 440735;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	840435	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (440987)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.567GB.  
  PID: 440987;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_noMT.bam` (441018,441019,441020,441021)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 6.567GB.  
  PID: 441020;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 441018;	Command: samtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 441019;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 441021;	Command: xargs;	Return code: 0;	Memory used: 0.064GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (441075)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 441075;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam` (441076)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.567GB.  
  PID: 441076;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


### Split BAM file (02-27 12:21:41) elapsed: 199.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam` (441100,441101)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:50. Running peak memory: 6.567GB.  
  PID: 441100;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 441101;	Command: samtools;	Return code: 0;	Memory used: 4.688GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE2.bam` (441301,441302)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:18. Running peak memory: 6.567GB.  
  PID: 441302;	Command: samtools;	Return code: 0;	Memory used: 3.629GB  
  PID: 441301;	Command: samtools;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam` (441716)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 6.567GB.  
  PID: 441716;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam` (441755,441756)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:43. Running peak memory: 6.567GB.  
  PID: 441755;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 441756;	Command: samtools;	Return code: 0;	Memory used: 4.448GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE2.bam` (442172,442173)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:09. Running peak memory: 6.567GB.  
  PID: 442172;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 442173;	Command: samtools;	Return code: 0;	Memory used: 3.566GB


### Calculate library complexity (02-27 12:33:00) elapsed: 679.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_out.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam` (442336)
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
Command completed. Elapsed time: 0:01:45. Running peak memory: 6.567GB.  
  PID: 442336;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam` (442436)
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
................................................._...................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:02:02. Running peak memory: 6.567GB.  
  PID: 442436;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_counts.txt` (442791)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 6.567GB.  
  PID: 442791;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_plot` (442835)
<pre>
Processing H9_PRO-seq_2
INFO: Found real counts for H9_PRO-seq_2 - Total (M): 19.4822 Unique (M): 18.471612

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.567GB.  
  PID: 442835;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `Library complexity`	QC_hg38/H9_PRO-seq_2_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_2_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.854	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 12:37:21) elapsed: 261.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam` (442854)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.567GB.  
  PID: 442854;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv` (442868)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/tmp_H9_PRO-seq_2_PE1_4koe83tc'
Processing with 8 cores...
Discarding 100 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270755v1']
Keeping 95 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 6.567GB.  
  PID: 442868;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.21GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	9741100.0	PEPPRO	_RES_

> `PBC2`	9741100.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_unmap.bam` (442917)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.567GB.  
  PID: 442917;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_temp.bam`

> `Unmapped_reads`	6664276	PEPPRO	_RES_

### Split BAM by strand (02-27 12:38:18) elapsed: 56.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam` (442957)
<pre>
</pre>
Command completed. Elapsed time: 0:01:13. Running peak memory: 6.567GB.  
  PID: 442957;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam` (443030)
<pre>
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 6.567GB.  
  PID: 443030;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 12:40:39) elapsed: 141.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (443315)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 443315;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt` (443316)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.567GB.  
  PID: 443316;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.497GB


> `TSS_coding_score`	33.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt` (443342)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.567GB.  
  PID: 443342;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.491GB


> `TSS_non-coding_score`	12.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt` (443376)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.567GB.  
  PID: 443376;	Command: Rscript;	Return code: 0;	Memory used: 0.125GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_2_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_2_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt` (443399,443400,443401,443402)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 443399;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 443401;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 443400;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 443402;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt` (443404)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 443404;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-27 12:41:05) elapsed: 26.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed` (443406,443407)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.567GB.  
  PID: 443406;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 443407;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed` (443411,443412)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 443411;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 443412;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSS_density.bed` (443414,443415,443416,443417)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.567GB.  
  PID: 443415;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 443417;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 443414;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 443416;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_gene_body_density.bed` (443446,443447,443448)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 6.567GB.  
  PID: 443448;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 443446;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB  
  PID: 443447;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed` (443486,443487,443488)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 443486;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 443488;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 443487;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	41.45	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed` (443496)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.567GB.  
  PID: 443496;	Command: Rscript;	Return code: 0;	Memory used: 0.262GB

> `Pause index`	QC_hg38/H9_PRO-seq_2_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_2_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_pause_index.bed` (443515)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 443515;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 12:42:16) elapsed: 71.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam`
18165727.0 7048811

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam`
18165727.0 6572286

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_gene_sort.bed` (443570,443571)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 443570;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 443571;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_gene_coverage.bed` (443575)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 6.567GB.  
  PID: 443575;	Command: bedtools;	Return code: 0;	Memory used: 0.066GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed.gz` (443609)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 443609;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed` (443610)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 443610;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 12:43:28) elapsed: 72.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/raw/hg38_annotations.bed` (443620)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 443620;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR"` (443622)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 443622;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed` (443623,443624,443625,443626)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 443623;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 443624;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 443626;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 443625;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed` (443629)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 443629;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_minus_coverage.bed` (443641)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 443641;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR"` (443659)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 443659;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed` (443660,443661,443662,443663)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 443660;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 443661;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 443663;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 443662;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed` (443666)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 443666;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_minus_coverage.bed` (443681)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 443681;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed` (443697,443698,443699,443700)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 443697;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 443698;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 443700;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 443699;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed` (443703)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 443703;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_minus_coverage.bed` (443716)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.567GB.  
  PID: 443716;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed` (443730,443731,443733,443734)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 6.567GB.  
  PID: 443730;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 443731;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 443734;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 443733;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed` (443741)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.567GB.  
  PID: 443741;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_minus_coverage.bed` (443949)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.567GB.  
  PID: 443949;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed` (443973,443974,443975,443976)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.567GB.  
  PID: 443973;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 443975;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 443974;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 443976;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed` (443979)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.567GB.  
  PID: 443979;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_minus_coverage.bed` (443994)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.567GB.  
  PID: 443994;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed` (444014,444015,444016,444017)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 444014;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 444015;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 444017;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 444016;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed` (444019)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.567GB.  
  PID: 444019;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_minus_coverage.bed` (444038)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.567GB.  
  PID: 444038;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region"` (444051)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 444051;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed` (444052,444053,444054,444055)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 444052;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 444054;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 444053;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 444055;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (444058)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.567GB.  
  PID: 444058;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed` (444072)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.567GB.  
  PID: 444072;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB


### Plot cFRiF/FRiF (02-27 12:46:55) elapsed: 207.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_2 -z 3099922541 -n 9944138 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (444104)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 6.567GB.  
  PID: 444104;	Command: Rscript;	Return code: 0;	Memory used: 0.442GB

> `cFRiF`	QC_hg38/H9_PRO-seq_2_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_2_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_PRO-seq_2 -z 3099922541 -n 9944138 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (444154)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.567GB.  
  PID: 444154;	Command: Rscript;	Return code: 0;	Memory used: 0.463GB

> `FRiF`	QC_hg38/H9_PRO-seq_2_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_2_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 12:48:05) elapsed: 70.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_exons_sort.bed` (444197,444198)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.567GB.  
  PID: 444198;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB  
  PID: 444197;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_introns_sort.bed` (444207,444208,444209)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.567GB.  
  PID: 444207;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 444209;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 444208;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_coverage.bed` (444216)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.567GB.  
  PID: 444216;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_coverage.bed` (444238)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 6.567GB.  
  PID: 444238;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.165727)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_rpkm.bed` (444276,444277,444278)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 444276;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 444278;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 444277;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.165727)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_rpkm.bed` (444281,444282,444283)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.567GB.  
  PID: 444281;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 444283;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 444282;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed` (444285,444286,444287)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 444285;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 444287;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 444286;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed --annotate` (444293)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.567GB.  
  PID: 444293;	Command: Rscript;	Return code: 0;	Memory used: 0.236GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_2_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_2_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/QC_hg38/H9_PRO-seq_2_exon_intron_ratios.bed` (444312)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.567GB.  
  PID: 444312;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-27 12:49:19) elapsed: 73.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam` (444320)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.567GB.  
  PID: 444320;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (444328)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_2_plus_cuttrace_5tr_64zz'
Processing with 2 cores...
stdin is empty of data
Discarding 112 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 83 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_exact_body_0-mer.bw'
Merging 83 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:40. Running peak memory: 6.567GB.  
  PID: 444328;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.169GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam` (446373)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.567GB.  
  PID: 446373;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (446386)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/aligned_hg38/H9_PRO-seq_2_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_2_minus_cuttrace_pxl7_jge'
Processing with 2 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr14_KI270723v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270755v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_exact_body_0-mer.bw'
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_PRO-seq_2/signal_hg38/H9_PRO-seq_2_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:54. Running peak memory: 6.567GB.  
  PID: 446386;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.508GB

Starting cleanup: 77 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:45:40
*  Total elapsed time (all runs):  7:15:29
*         Peak memory (this run):  6.5669 GB
*        Pipeline completed time: 2020-02-27 13:11:12
