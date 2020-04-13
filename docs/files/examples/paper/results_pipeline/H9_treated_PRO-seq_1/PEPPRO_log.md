### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name H9_treated_PRO-seq_1 --genome hg38 --input /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz --protocol PRO --umi-len 8 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj38-14c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/
*  Pipeline started at:   (02-27 09:26:04) elapsed: 1.0 _TIME_

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
*              `input`:  `['/project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz']`
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

### Merge/link and fastq conversion:  (02-27 09:26:05) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz` (333749)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 333749;	Command: ln;	Return code: 0;	Memory used: 0.002GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz` (333750)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.002GB.  
  PID: 333750;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1.fastq`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1.fastq` (333751)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 0.003GB.  
  PID: 333751;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2.fastq` (136808)
<pre>
</pre>
Command completed. Elapsed time: 0:00:57. Running peak memory: 0.003GB.  
  PID: 136808;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	111300168	PEPPRO	_RES_

> `Fastq_reads`	111300168	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz']

### FASTQ processing:  (02-27 09:29:37) elapsed: 212.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt` (407832)
<pre>
</pre>
Command completed. Elapsed time: 0:01:19. Running peak memory: 0.279GB.  
  PID: 407832;	Command: cutadapt;	Return code: 0;	Memory used: 0.279GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (408178,408179)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 0.279GB.  
  PID: 408178;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 408179;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	29119715	PEPPRO	_RES_

> `Trim_loss_rate`	73.84	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (408216)
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
Command completed. Elapsed time: 0:01:17. Running peak memory: 0.279GB.  
  PID: 408216;	Command: fastqc;	Return code: 0;	Memory used: 0.174GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_dedup.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq` (408315)
<pre>
[INFO][0m 2891628 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 2.527GB.  
  PID: 408315;	Command: seqkit;	Return code: 0;	Memory used: 2.527GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq` (408370,408371)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 2.527GB.  
  PID: 408370;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 408371;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	43631018.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	25201427.0	PEPPRO	_RES_

> `Duplicate_reads`	2891628.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	45.2855	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/processed_R1.flag` (408452)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.527GB.  
  PID: 408452;	Command: touch;	Return code: 0;	Memory used: 0.002GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R2_cutadapt.txt` (408455)
<pre>
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 2.527GB.  
  PID: 408455;	Command: cutadapt;	Return code: 0;	Memory used: 0.264GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (408762,408763)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 2.527GB.  
  PID: 408762;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 408763;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	58239430	PEPPRO	_RES_

> `Trim_loss_rate`	47.67	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (408807)
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
Command completed. Elapsed time: 0:01:19. Running peak memory: 2.527GB.  
  PID: 408807;	Command: fastqc;	Return code: 0;	Memory used: 0.175GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (408933)
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
Command completed. Elapsed time: 0:01:14. Running peak memory: 2.527GB.  
  PID: 408933;	Command: fastqc;	Return code: 0;	Memory used: 0.168GB

> `FastQC report r2`	fastqc/H9_treated_PRO-seq_1_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.histogram`  

> `fastq_pair -t 100170151 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq` (409030)
<pre>
Left paired: 30307647		Right paired: 30307647
Left single: 141010		Right single: 2158148
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:44. Running peak memory: 6.179GB.  
  PID: 409030;	Command: fastq_pair;	Return code: 0;	Memory used: 6.179GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq.paired.fq -o H9_treated_PRO-seq_1 -d /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt` (409421)
<pre>
</pre>
Command completed. Elapsed time: 0:01:24. Running peak memory: 6.179GB.  
  PID: 409421;	Command: flash;	Return code: 0;	Memory used: 0.072GB


### Plot adapter insertion distribution (02-27 09:43:01) elapsed: 804.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt -u 8` (409651)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.179GB.  
  PID: 409651;	Command: Rscript;	Return code: 0;	Memory used: 0.107GB

> `Adapter insertion distribution`	cutadapt/H9_treated_PRO-seq_1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_treated_PRO-seq_1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:43:09) elapsed: 8.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `Degradation_ratio`	1.158	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq` (409681)
<pre>
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 6.179GB.  
  PID: 409681;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/processed_R2.flag` (409700)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.179GB.  
  PID: 409700;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/repaired.flag`  

> `fastq_pair -t 100170151 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (409701)
<pre>
Left paired: 28948010		Right paired: 28948010
Left single: 171705		Right single: 2245320
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:10. Running peak memory: 6.48GB.  
  PID: 409701;	Command: fastq_pair;	Return code: 0;	Memory used: 6.48GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (410044)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 410044;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (410047)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 410047;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/repaired.flag` (410049)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 410049;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/dups_repaired.flag`  

> `fastq_pair -t 100170151 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq` (410050)
<pre>
Left paired: 26293117		Right paired: 26293117
Left single: 111289		Right single: 4900213
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:57. Running peak memory: 6.48GB.  
  PID: 410050;	Command: fastq_pair;	Return code: 0;	Memory used: 6.176GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq` (410167)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 410167;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq` (410168)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 410168;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/dups_repaired.flag` (410170)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 410170;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-27 09:47:37) elapsed: 268.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:47:37) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2` (410171)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 410171;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq` (410172)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
28948010 reads; of these:
  28948010 (100.00%) were unpaired; of these:
    25898407 (89.47%) aligned 0 times
    3049603 (10.53%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.53% overall alignment rate

> `Aligned_reads_human_rDNA`	6099206.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.47	PEPPRO	_RES_

### Map to human_rDNA (02-27 09:51:09) elapsed: 211.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2` (410631)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 410631;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R2.fq` (410632)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
3049603 reads skipped
0 reads lost

### Map to genome (02-27 09:54:22) elapsed: 193.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_1 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/tmp_9nne37s -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam` (410837,410838,410839)
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
Command completed. Elapsed time: 1:04:33. Running peak memory: 6.48GB.  
  PID: 410838;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 410837;	Command: bowtie2;	Return code: 0;	Memory used: 3.687GB  
  PID: 410839;	Command: samtools;	Return code: 0;	Memory used: 0.872GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (417738)
<pre>
</pre>
Command completed. Elapsed time: 0:01:51. Running peak memory: 6.48GB.  
  PID: 417738;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	44081295	PEPPRO	_RES_

> `QC_filtered_reads`	26072183	PEPPRO	_RES_

> `Aligned_reads`	18009111.5	PEPPRO	_RES_

> `Alignment_rate`	30.92	PEPPRO	_RES_

> `Total_efficiency`	16.18	PEPPRO	_RES_

> `Read_depth`	3.41	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_1 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/tmp_9nne37s -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam` (419333,419339,419340)
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
Command completed. Elapsed time: 0:57:02. Running peak memory: 6.48GB.  
  PID: 419333;	Command: bowtie2;	Return code: 0;	Memory used: 3.663GB  
  PID: 419339;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 419340;	Command: samtools;	Return code: 0;	Memory used: 0.872GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam` (425327)
<pre>
</pre>
Command completed. Elapsed time: 0:01:35. Running peak memory: 6.48GB.  
  PID: 425327;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


### Compress all unmapped read files (02-27 12:11:31) elapsed: 8229.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq` (425670)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 6.48GB.  
  PID: 425670;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq` (425718)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.48GB.  
  PID: 425718;	Command: pigz;	Return code: 0;	Memory used: 0.006GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam` (425755)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 6.48GB.  
  PID: 425755;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	655778	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (425803)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.48GB.  
  PID: 425803;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_noMT.bam` (425830,425831,425832,425833)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 6.48GB.  
  PID: 425832;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 425830;	Command: samtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 425831;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 425833;	Command: xargs;	Return code: 0;	Memory used: 0.068GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (425882)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 425882;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (425883)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.48GB.  
  PID: 425883;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


### Split BAM file (02-27 12:14:44) elapsed: 193.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam` (425912,425913)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:40. Running peak memory: 6.48GB.  
  PID: 425912;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 425913;	Command: samtools;	Return code: 0;	Memory used: 4.884GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE2.bam` (426311,426312)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:10. Running peak memory: 6.48GB.  
  PID: 426311;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 426312;	Command: samtools;	Return code: 0;	Memory used: 3.784GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam` (426798)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 6.48GB.  
  PID: 426798;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam` (426840,426841)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:32. Running peak memory: 6.48GB.  
  PID: 426840;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 426841;	Command: samtools;	Return code: 0;	Memory used: 4.622GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE2.bam` (427021,427022)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:02. Running peak memory: 6.48GB.  
  PID: 427021;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 427022;	Command: samtools;	Return code: 0;	Memory used: 3.576GB


### Calculate library complexity (02-27 12:25:25) elapsed: 641.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_out.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam` (427385)
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
Command completed. Elapsed time: 0:01:46. Running peak memory: 6.48GB.  
  PID: 427385;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam` (427494)
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
.............................................._........._.............................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:53. Running peak memory: 6.48GB.  
  PID: 427494;	Command: preseq;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_counts.txt` (427728)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.48GB.  
  PID: 427728;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_plot` (427760)
<pre>
Processing H9_treated_PRO-seq_1
INFO: Found real counts for H9_treated_PRO-seq_1 - Total (M): 19.73867 Unique (M): 18.664481

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.48GB.  
  PID: 427760;	Command: Rscript;	Return code: 0;	Memory used: 0.221GB

> `Library complexity`	QC_hg38/H9_treated_PRO-seq_1_preseq_plot.pdf	Library complexity	QC_hg38/H9_treated_PRO-seq_1_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8095	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 12:29:37) elapsed: 253.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam` (427779)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.48GB.  
  PID: 427779;	Command: samtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv` (427800)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/tmp_H9_treated_PRO-seq_1_PE1_xv9yaelj'
Processing with 8 cores...
Discarding 94 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000214v1']
Keeping 101 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.48GB.  
  PID: 427800;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.895GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	9869335.0	PEPPRO	_RES_

> `PBC2`	9869335.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_unmap.bam` (428089)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 6.48GB.  
  PID: 428089;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam`

> `Unmapped_reads`	7715519	PEPPRO	_RES_

### Split BAM by strand (02-27 12:30:37) elapsed: 59.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam` (428126)
<pre>
</pre>
Command completed. Elapsed time: 0:01:10. Running peak memory: 6.48GB.  
  PID: 428126;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam` (428205)
<pre>
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 6.48GB.  
  PID: 428205;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 12:32:56) elapsed: 139.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (428273)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 428273;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt` (428274)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 6.48GB.  
  PID: 428274;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.488GB


> `TSS_coding_score`	53.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt` (428302)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.48GB.  
  PID: 428302;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.484GB


> `TSS_non-coding_score`	17.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt` (428331)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.48GB.  
  PID: 428331;	Command: Rscript;	Return code: 0;	Memory used: 0.264GB

> `TSS enrichment`	QC_hg38/H9_treated_PRO-seq_1_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_treated_PRO-seq_1_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt` (428350,428351,428352,428353)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 428350;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 428352;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 428351;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 428353;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt` (428355)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 428355;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (02-27 12:33:19) elapsed: 24.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed` (428357,428358)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.48GB.  
  PID: 428357;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 428358;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed` (428361,428362)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 428361;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 428362;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSS_density.bed` (428365,428366,428367,428368)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.48GB.  
  PID: 428366;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 428368;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 428365;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 428367;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_gene_body_density.bed` (428401,428402,428403)
<pre>
</pre>
Command completed. Elapsed time: 0:00:43. Running peak memory: 6.48GB.  
  PID: 428401;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB  
  PID: 428403;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 428402;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed` (428443,428445,428446)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 428443;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 428446;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 428445;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	75.18	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed` (428451)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.48GB.  
  PID: 428451;	Command: Rscript;	Return code: 0;	Memory used: 0.253GB

> `Pause index`	QC_hg38/H9_treated_PRO-seq_1_pause_index.pdf	Pause index	QC_hg38/H9_treated_PRO-seq_1_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed` (428472)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 428472;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 12:34:40) elapsed: 80.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam`
18009111.5 7282779

> `Plus_FRiP`	0.4	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam`
18009111.5 6888359

> `Minus_FRiP`	0.38	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_gene_sort.bed` (428741,428742)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 428741;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 428742;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_gene_coverage.bed` (428745)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 6.48GB.  
  PID: 428745;	Command: bedtools;	Return code: 0;	Memory used: 0.06GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed.gz` (428788)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 428788;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed` (428789)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 428789;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 12:35:57) elapsed: 78.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed` (428798)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 428798;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR"` (428800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 428800;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed` (428801,428802,428803,428804)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 428801;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 428802;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 428804;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 428803;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed` (428806)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.48GB.  
  PID: 428806;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_minus_coverage.bed` (428822)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.48GB.  
  PID: 428822;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR"` (428867)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 428867;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed` (428868,428869,428870,428871)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 428868;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 428869;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 428871;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 428870;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed` (428873)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.48GB.  
  PID: 428873;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_minus_coverage.bed` (428889)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.48GB.  
  PID: 428889;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed` (428908,428909,428910,428911)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 428908;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 428909;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 428911;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB  
  PID: 428910;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed` (428913)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.48GB.  
  PID: 428913;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_minus_coverage.bed` (428926)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.48GB.  
  PID: 428926;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed` (428941,428942,428943,428944)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.48GB.  
  PID: 428941;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 428942;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 428944;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 428943;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed` (428949)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.48GB.  
  PID: 428949;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_minus_coverage.bed` (428967)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.48GB.  
  PID: 428967;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed` (428986,428987,428988,428989)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 428986;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 428988;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 428987;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 428989;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed` (428992)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.48GB.  
  PID: 428992;	Command: bedtools;	Return code: 0;	Memory used: 0.076GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_minus_coverage.bed` (429011)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 6.48GB.  
  PID: 429011;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed` (429027,429028,429029,429030)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 429027;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 429028;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 429030;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 429029;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed` (429032)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.48GB.  
  PID: 429032;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_minus_coverage.bed` (429053)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.48GB.  
  PID: 429053;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region"` (429069)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 429069;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed` (429070,429071,429072,429073)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 429070;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 429072;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 429071;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 429073;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (429076)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.48GB.  
  PID: 429076;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed` (429092)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.48GB.  
  PID: 429092;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


### Plot cFRiF/FRiF (02-27 12:39:28) elapsed: 211.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_1 -z 3099922541 -n 10005525 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (429117)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 6.48GB.  
  PID: 429117;	Command: Rscript;	Return code: 0;	Memory used: 0.445GB

> `cFRiF`	QC_hg38/H9_treated_PRO-seq_1_cFRiF.pdf	cFRiF	QC_hg38/H9_treated_PRO-seq_1_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_1 -z 3099922541 -n 10005525 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (429391)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 6.48GB.  
  PID: 429391;	Command: Rscript;	Return code: 0;	Memory used: 0.476GB

> `FRiF`	QC_hg38/H9_treated_PRO-seq_1_FRiF.pdf	FRiF	QC_hg38/H9_treated_PRO-seq_1_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 12:40:37) elapsed: 68.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_exons_sort.bed` (429434,429435)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.48GB.  
  PID: 429435;	Command: bedtools;	Return code: 0;	Memory used: 0.085GB  
  PID: 429434;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_introns_sort.bed` (429445,429446,429447)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.48GB.  
  PID: 429445;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 429447;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 429446;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_coverage.bed` (429458)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.48GB.  
  PID: 429458;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_coverage.bed` (429485)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 6.48GB.  
  PID: 429485;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.0091115)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_rpkm.bed` (429523,429524,429525)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 429523;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 429525;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 429524;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.0091115)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_rpkm.bed` (429528,429529,429530)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.48GB.  
  PID: 429528;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 429530;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 429529;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed` (429532,429533,429534)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 429532;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 429534;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 429533;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.19	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed --annotate` (429540)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.48GB.  
  PID: 429540;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `mRNA contamination`	QC_hg38/H9_treated_PRO-seq_1_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_treated_PRO-seq_1_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed` (429560)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.48GB.  
  PID: 429560;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (02-27 12:41:55) elapsed: 79.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam` (429568)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.48GB.  
  PID: 429568;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (429575)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_1_plus_cuttrace_bsfi2yam'
Processing with 2 cores...
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270752v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:19. Running peak memory: 6.48GB.  
  PID: 429575;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.549GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam` (431435)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.48GB.  
  PID: 431435;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (431443)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_1_minus_cuttrace_hm1p2zhr'
Processing with 2 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 106 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr17_KI270730v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1']
Keeping 89 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 89 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_exact_body_0-mer.bw'
Merging 89 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:11:01. Running peak memory: 6.48GB.  
  PID: 431443;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.549GB

Starting cleanup: 77 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:37:30
*  Total elapsed time (all runs):  7:09:26
*         Peak memory (this run):  6.4799 GB
*        Pipeline completed time: 2020-02-27 13:03:34
