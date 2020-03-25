### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name H9_treated_PRO-seq_2 --genome hg38 --input /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep2_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep2_PE2.fastq.gz --protocol PRO --umi-len 8 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj38-14c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/
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
*              `input`:  `['/project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep2_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep2_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_treated_PRO-seq_2`
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

Local input file: /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep2_PE1.fastq.gz
Local input file: /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep2_PE2.fastq.gz

> `File_mb`	2737.32	PEPPRO	_RES_

> `Read_type`	paired	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-27 09:26:04) elapsed: 0.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep2_PE1.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R1.fastq.gz` (420900)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 420900;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/H9_200nM_romidepsin_rep2_PE2.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R2.fastq.gz` (420901)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 420901;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1.fastq`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R1.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1.fastq` (420902)
<pre>
</pre>
Command completed. Elapsed time: 0:00:55. Running peak memory: 0.003GB.  
  PID: 420902;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R2.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2.fastq` (227640)
<pre>
</pre>
Command completed. Elapsed time: 0:00:58. Running peak memory: 0.003GB.  
  PID: 227640;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	114083518	PEPPRO	_RES_

> `Fastq_reads`	114083518	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R2.fastq.gz']

### FASTQ processing:  (02-27 09:29:47) elapsed: 223.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R1_cutadapt.txt` (56051)
<pre>
</pre>
Command completed. Elapsed time: 0:01:16. Running peak memory: 0.285GB.  
  PID: 56051;	Command: cutadapt;	Return code: 0;	Memory used: 0.285GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq` (56384,56385)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 0.285GB.  
  PID: 56384;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 56385;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	28674037	PEPPRO	_RES_

> `Trim_loss_rate`	74.87	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq` (56417)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_treated_PRO-seq_2_R1_processed.fastq
Approx 5% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 10% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 15% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 20% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 25% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 30% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 35% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 40% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 45% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 50% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 55% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 60% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 65% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 70% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 75% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 80% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 85% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 90% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 95% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 100% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Analysis complete for H9_treated_PRO-seq_2_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 0.285GB.  
  PID: 56417;	Command: fastqc;	Return code: 0;	Memory used: 0.18GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_2_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_dedup.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq` (56502)
<pre>
[INFO][0m 2388789 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 1.968GB.  
  PID: 56502;	Command: seqkit;	Return code: 0;	Memory used: 1.968GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq` (56562,56563)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 1.968GB.  
  PID: 56562;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 56563;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	43984199.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	27614116.0	PEPPRO	_RES_

> `Duplicate_reads`	2388789.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	48.4104	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/processed_R1.flag` (56631)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.968GB.  
  PID: 56631;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R2_cutadapt.txt` (56634)
<pre>
</pre>
Command completed. Elapsed time: 0:01:16. Running peak memory: 1.968GB.  
  PID: 56634;	Command: cutadapt;	Return code: 0;	Memory used: 0.271GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq` (56927,56929)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 1.968GB.  
  PID: 56927;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 56929;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	57348074	PEPPRO	_RES_

> `Trim_loss_rate`	49.73	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq` (56975)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_treated_PRO-seq_2_R1_processed.fastq
Approx 5% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 10% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 15% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 20% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 25% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 30% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 35% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 40% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 45% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 50% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 55% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 60% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 65% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 70% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 75% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 80% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 85% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 90% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 95% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Approx 100% complete for H9_treated_PRO-seq_2_R1_processed.fastq
Analysis complete for H9_treated_PRO-seq_2_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 1.968GB.  
  PID: 56975;	Command: fastqc;	Return code: 0;	Memory used: 0.178GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_2_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq` (57051)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 5% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 10% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 15% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 20% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 25% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 30% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 35% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 40% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 45% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 50% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 55% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 60% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 65% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 70% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 75% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 80% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 85% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 90% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Approx 95% complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
Analysis complete for H9_treated_PRO-seq_2_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 1.968GB.  
  PID: 57051;	Command: fastqc;	Return code: 0;	Memory used: 0.172GB

> `FastQC report r2`	fastqc/H9_treated_PRO-seq_2_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.histogram`  

> `fastq_pair -t 102675166 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq` (57128)
<pre>
Left paired: 29252476		Right paired: 29252476
Left single: 175167		Right single: 2238349
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:53. Running peak memory: 6.557GB.  
  PID: 57128;	Command: fastq_pair;	Return code: 0;	Memory used: 6.557GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq.paired.fq -o H9_treated_PRO-seq_2 -d /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt` (57547)
<pre>
</pre>
Command completed. Elapsed time: 0:01:18. Running peak memory: 6.557GB.  
  PID: 57547;	Command: flash;	Return code: 0;	Memory used: 0.127GB


### Plot adapter insertion distribution (02-27 09:42:51) elapsed: 784.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt -u 8` (57767)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.557GB.  
  PID: 57767;	Command: Rscript;	Return code: 0;	Memory used: 0.107GB

> `Adapter insertion distribution`	cutadapt/H9_treated_PRO-seq_2_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_treated_PRO-seq_2_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:42:59) elapsed: 8.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`

> `Degradation_ratio`	1.0037	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq` (57799)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.557GB.  
  PID: 57799;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/processed_R2.flag` (57822)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 57822;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/repaired.flag`  

> `fastq_pair -t 102675166 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq` (57823)
<pre>
Left paired: 28465823		Right paired: 28465823
Left single: 208214		Right single: 2295865
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:40. Running peak memory: 6.557GB.  
  PID: 57823;	Command: fastq_pair;	Return code: 0;	Memory used: 6.004GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq` (58177)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 58177;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq` (58179)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 58179;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/repaired.flag` (58180)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 58180;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/dups_repaired.flag`  

> `fastq_pair -t 102675166 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq` (58181)
<pre>
Left paired: 26273211		Right paired: 26273211
Left single: 141535		Right single: 4488477
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:42. Running peak memory: 6.557GB.  
  PID: 58181;	Command: fastq_pair;	Return code: 0;	Memory used: 4.33GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq` (58272)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 58272;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq` (58273)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 58273;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/dups_repaired.flag` (58274)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 58274;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-27 09:47:50) elapsed: 291.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:47:50) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_bt2` (58275)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 58275;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq` (58276)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_2 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
28465823 reads; of these:
  28465823 (100.00%) were unpaired; of these:
    25736513 (90.41%) aligned 0 times
    2729310 (9.59%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.59% overall alignment rate

> `Aligned_reads_human_rDNA`	5458620.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.52	PEPPRO	_RES_

### Map to human_rDNA (02-27 09:51:13) elapsed: 203.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_dups_bt2` (58698)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 58698;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_dups_R2.fq` (58699)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_2 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2729310 reads skipped
0 reads lost

### Map to genome (02-27 09:54:36) elapsed: 202.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_2 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/tmp4i5ti6yr -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam` (58897,58898,58899)
<pre>
2306063 reads skipped
0 reads lost
25736513 reads; of these:
  25736513 (100.00%) were paired; of these:
    10470203 (40.68%) aligned concordantly 0 times
    12826144 (49.84%) aligned concordantly exactly 1 time
    2440166 (9.48%) aligned concordantly >1 times
    ----
    10470203 pairs aligned concordantly 0 times; of these:
      2441340 (23.32%) aligned discordantly 1 time
    ----
    8028863 pairs aligned 0 times concordantly or discordantly; of these:
      16057726 mates make up the pairs; of these:
        6731776 (41.92%) aligned 0 times
        3441445 (21.43%) aligned exactly 1 time
        5884505 (36.65%) aligned >1 times
86.92% overall alignment rate
[bam_sort_core] merging from 14 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:59:44. Running peak memory: 6.557GB.  
  PID: 58897;	Command: bowtie2;	Return code: 0;	Memory used: 3.69GB  
  PID: 58898;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 58899;	Command: samtools;	Return code: 0;	Memory used: 0.872GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam` (64877)
<pre>
</pre>
Command completed. Elapsed time: 0:01:31. Running peak memory: 6.557GB.  
  PID: 64877;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	44741250	PEPPRO	_RES_

> `QC_filtered_reads`	25860518	PEPPRO	_RES_

> `Aligned_reads`	18880731.5	PEPPRO	_RES_

> `Alignment_rate`	32.92	PEPPRO	_RES_

> `Total_efficiency`	16.55	PEPPRO	_RES_

> `Read_depth`	3.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_2 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/tmp4i5ti6yr -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp_dups.bam` (66318,66325,66326)
<pre>
23967148 reads; of these:
  23967148 (100.00%) were paired; of these:
    9441036 (39.39%) aligned concordantly 0 times
    12234981 (51.05%) aligned concordantly exactly 1 time
    2291131 (9.56%) aligned concordantly >1 times
    ----
    9441036 pairs aligned concordantly 0 times; of these:
      2335210 (24.73%) aligned discordantly 1 time
    ----
    7105826 pairs aligned 0 times concordantly or discordantly; of these:
      14211652 mates make up the pairs; of these:
        6043583 (42.53%) aligned 0 times
        3281410 (23.09%) aligned exactly 1 time
        4886659 (34.38%) aligned >1 times
87.39% overall alignment rate
[bam_sort_core] merging from 13 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:53:01. Running peak memory: 6.557GB.  
  PID: 66325;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 66318;	Command: bowtie2;	Return code: 0;	Memory used: 3.689GB  
  PID: 66326;	Command: samtools;	Return code: 0;	Memory used: 0.879GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp_dups.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam` (71719)
<pre>
</pre>
Command completed. Elapsed time: 0:01:22. Running peak memory: 6.557GB.  
  PID: 71719;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (02-27 12:02:40) elapsed: 7685.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq` (71805)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.557GB.  
  PID: 71805;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R1.fq` (71853)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.557GB.  
  PID: 71853;	Command: pigz;	Return code: 0;	Memory used: 0.01GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam` (71891)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 6.557GB.  
  PID: 71891;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	670904	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam` (71928)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.557GB.  
  PID: 71928;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_noMT.bam` (71957,71958,71959,71960)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 6.557GB.  
  PID: 71959;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 71957;	Command: samtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 71958;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 71960;	Command: xargs;	Return code: 0;	Memory used: 0.06GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam` (72215)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 72215;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam` (72216)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.557GB.  
  PID: 72216;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


### Split BAM file (02-27 12:05:54) elapsed: 194.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam` (72243,72244)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:39. Running peak memory: 6.557GB.  
  PID: 72243;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 72244;	Command: samtools;	Return code: 0;	Memory used: 5.036GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE2.bam` (72409,72410)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:15. Running peak memory: 6.557GB.  
  PID: 72409;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 72410;	Command: samtools;	Return code: 0;	Memory used: 4.082GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp_dups.bam` (72865)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 6.557GB.  
  PID: 72865;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE1.bam` (72902,72903)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:32. Running peak memory: 6.557GB.  
  PID: 72902;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 72903;	Command: samtools;	Return code: 0;	Memory used: 4.879GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE2.bam` (73073,73074)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:11. Running peak memory: 6.557GB.  
  PID: 73073;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 73074;	Command: samtools;	Return code: 0;	Memory used: 3.946GB


### Calculate library complexity (02-27 12:16:48) elapsed: 654.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_out.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE1.bam` (73432)
<pre>
BAM_INPUT
TOTAL READS     = 19595024
COUNTS_SUM      = 19595024
DISTINCT READS  = 1.48576e+07
DISTINCT COUNTS = 296
MAX COUNT       = 22664
COUNTS OF 1     = 1.29529e+07
OBSERVED COUNTS (22665)
1	12952887
2	1171763
3	313507
4	140501
5	77968
6	49441
7	32956
8	23332
9	17067
10	12863
11	9993
12	8220
13	6438
14	5374
15	4272
16	3535
17	3085
18	2617
19	2187
20	1897
21	1621
22	1418
23	1206
24	1109
25	972
26	869
27	788
28	680
29	637
30	560
31	514
32	459
33	411
34	377
35	353
36	323
37	323
38	279
39	249
40	227
41	222
42	223
43	172
44	153
45	195
46	137
47	127
48	127
49	119
50	121
51	114
52	120
53	100
54	84
55	90
56	84
57	60
58	80
59	77
60	77
61	68
62	80
63	66
64	44
65	54
66	49
67	32
68	53
69	48
70	41
71	32
72	39
73	32
74	29
75	29
76	27
77	31
78	22
79	29
80	27
81	31
82	32
83	33
84	21
85	27
86	21
87	8
88	19
89	15
90	15
91	12
92	20
93	18
94	14
95	14
96	21
97	8
98	16
99	9
100	15
101	8
102	13
103	19
104	13
105	14
106	6
107	11
108	12
109	13
110	10
111	8
112	13
113	7
114	9
115	3
116	11
117	9
118	10
119	5
120	5
121	7
122	4
123	9
124	6
125	9
126	4
127	6
128	7
129	6
130	8
131	3
132	7
133	8
134	7
135	3
136	7
137	5
138	8
139	5
140	10
141	4
142	3
143	3
144	4
145	2
146	2
147	7
148	4
149	4
150	3
151	3
153	4
154	2
155	6
156	2
157	3
158	4
160	5
161	2
162	2
163	1
164	1
165	1
166	6
167	3
168	1
169	5
170	1
171	5
172	4
173	1
174	4
175	5
176	1
177	4
178	1
179	4
180	2
181	4
183	1
184	4
185	6
188	1
189	2
190	1
191	1
194	1
195	2
196	3
197	3
198	2
199	1
200	1
202	2
203	4
206	2
207	1
208	2
209	1
211	1
212	1
213	3
214	3
216	2
217	1
221	1
222	2
224	1
225	1
227	1
228	1
231	2
234	4
236	1
237	1
240	1
242	2
246	1
247	1
250	2
253	1
254	1
255	1
256	1
257	1
259	1
260	1
263	3
265	1
266	1
267	1
268	2
273	1
274	1
277	1
279	1
280	1
282	2
283	1
286	1
287	1
288	1
290	1
292	1
295	3
297	1
301	1
305	1
309	1
321	1
334	1
338	1
343	1
346	1
360	1
364	1
366	1
380	1
382	1
388	1
389	1
398	1
401	1
409	1
410	1
433	1
450	1
484	1
485	1
492	1
498	1
501	1
502	1
558	1
607	1
609	1
643	1
718	1
762	1
781	1
850	1
863	1
932	1
1085	1
1090	1
1227	1
1252	1
1420	1
1536	1
1567	1
1836	1
5591	1
5743	1
7315	1
15219	1
22664	1

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
Command completed. Elapsed time: 0:01:49. Running peak memory: 6.557GB.  
  PID: 73432;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE1.bam` (73525)
<pre>
BAM_INPUT
TOTAL READS     = 19595024
DISTINCT READS  = 1.48576e+07
DISTINCT COUNTS = 296
MAX COUNT       = 22664
COUNTS OF 1     = 1.29529e+07
MAX TERMS       = 100
OBSERVED COUNTS (22665)
1	12952887
2	1171763
3	313507
4	140501
5	77968
6	49441
7	32956
8	23332
9	17067
10	12863
11	9993
12	8220
13	6438
14	5374
15	4272
16	3535
17	3085
18	2617
19	2187
20	1897
21	1621
22	1418
23	1206
24	1109
25	972
26	869
27	788
28	680
29	637
30	560
31	514
32	459
33	411
34	377
35	353
36	323
37	323
38	279
39	249
40	227
41	222
42	223
43	172
44	153
45	195
46	137
47	127
48	127
49	119
50	121
51	114
52	120
53	100
54	84
55	90
56	84
57	60
58	80
59	77
60	77
61	68
62	80
63	66
64	44
65	54
66	49
67	32
68	53
69	48
70	41
71	32
72	39
73	32
74	29
75	29
76	27
77	31
78	22
79	29
80	27
81	31
82	32
83	33
84	21
85	27
86	21
87	8
88	19
89	15
90	15
91	12
92	20
93	18
94	14
95	14
96	21
97	8
98	16
99	9
100	15
101	8
102	13
103	19
104	13
105	14
106	6
107	11
108	12
109	13
110	10
111	8
112	13
113	7
114	9
115	3
116	11
117	9
118	10
119	5
120	5
121	7
122	4
123	9
124	6
125	9
126	4
127	6
128	7
129	6
130	8
131	3
132	7
133	8
134	7
135	3
136	7
137	5
138	8
139	5
140	10
141	4
142	3
143	3
144	4
145	2
146	2
147	7
148	4
149	4
150	3
151	3
153	4
154	2
155	6
156	2
157	3
158	4
160	5
161	2
162	2
163	1
164	1
165	1
166	6
167	3
168	1
169	5
170	1
171	5
172	4
173	1
174	4
175	5
176	1
177	4
178	1
179	4
180	2
181	4
183	1
184	4
185	6
188	1
189	2
190	1
191	1
194	1
195	2
196	3
197	3
198	2
199	1
200	1
202	2
203	4
206	2
207	1
208	2
209	1
211	1
212	1
213	3
214	3
216	2
217	1
221	1
222	2
224	1
225	1
227	1
228	1
231	2
234	4
236	1
237	1
240	1
242	2
246	1
247	1
250	2
253	1
254	1
255	1
256	1
257	1
259	1
260	1
263	3
265	1
266	1
267	1
268	2
273	1
274	1
277	1
279	1
280	1
282	2
283	1
286	1
287	1
288	1
290	1
292	1
295	3
297	1
301	1
305	1
309	1
321	1
334	1
338	1
343	1
346	1
360	1
364	1
366	1
380	1
382	1
388	1
389	1
398	1
401	1
409	1
410	1
433	1
450	1
484	1
485	1
492	1
498	1
501	1
502	1
558	1
607	1
609	1
643	1
718	1
762	1
781	1
850	1
863	1
932	1
1085	1
1090	1
1227	1
1252	1
1420	1
1536	1
1567	1
1836	1
5591	1
5743	1
7315	1
15219	1
22664	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
.........._....._........................_..........._..................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:58. Running peak memory: 6.557GB.  
  PID: 73525;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_counts.txt` (73857)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.557GB.  
  PID: 73857;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_plot` (73889)
<pre>
Processing H9_treated_PRO-seq_2
INFO: Found real counts for H9_treated_PRO-seq_2 - Total (M): 20.247016 Unique (M): 19.595024

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.557GB.  
  PID: 73889;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `Library complexity`	QC_hg38/H9_treated_PRO-seq_2_preseq_plot.pdf	Library complexity	QC_hg38/H9_treated_PRO-seq_2_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8186	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 12:21:09) elapsed: 261.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam` (73909)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.557GB.  
  PID: 73909;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_bamQC.tsv` (73922)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/tmp_H9_treated_PRO-seq_2_PE1_ficf2_sf'
Processing with 8 cores...
Discarding 98 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270757v1']
Keeping 97 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270372v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 6.557GB.  
  PID: 73922;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.06GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	10123508.0	PEPPRO	_RES_

> `PBC2`	10123508.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_unmap.bam` (73967)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.557GB.  
  PID: 73967;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam`

> `Unmapped_reads`	6731776	PEPPRO	_RES_

### Split BAM by strand (02-27 12:22:06) elapsed: 57.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam` (74004)
<pre>
</pre>
Command completed. Elapsed time: 0:01:11. Running peak memory: 6.557GB.  
  PID: 74004;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam` (74080)
<pre>
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 6.557GB.  
  PID: 74080;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 12:24:26) elapsed: 140.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (74141)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74141;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_plus_TssEnrichment.txt` (74142)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.557GB.  
  PID: 74142;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.611GB


> `TSS_coding_score`	57.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_minus_TssEnrichment.txt` (74168)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.557GB.  
  PID: 74168;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.612GB


> `TSS_non-coding_score`	16.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_minus_TssEnrichment.txt` (74193)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 6.557GB.  
  PID: 74193;	Command: Rscript;	Return code: 0;	Memory used: 0.343GB

> `TSS enrichment`	QC_hg38/H9_treated_PRO-seq_2_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_treated_PRO-seq_2_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt` (74215,74216,74217,74218)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74215;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 74217;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 74216;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 74218;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt` (74220)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74220;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-27 12:24:48) elapsed: 22.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed` (74222,74223)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.557GB.  
  PID: 74222;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 74223;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed` (74227,74228)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74227;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 74228;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_TSS_density.bed` (74230,74231,74232,74233)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.557GB.  
  PID: 74231;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 74233;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 74230;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 74232;	Command: sort;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_gene_body_density.bed` (74470,74471,74472)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.557GB.  
  PID: 74472;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 74470;	Command: bedtools;	Return code: 0;	Memory used: 0.058GB  
  PID: 74471;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed` (74503,74504,74505)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74503;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 74505;	Command: env;	Return code: 0;	Memory used: 0.007GB  
  PID: 74504;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	70.88	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed` (74510)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.557GB.  
  PID: 74510;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Pause index`	QC_hg38/H9_treated_PRO-seq_2_pause_index.pdf	Pause index	QC_hg38/H9_treated_PRO-seq_2_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed` (74530)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74530;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 12:25:57) elapsed: 69.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam`
18880731.5 7474169

> `Plus_FRiP`	0.4	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam`
18880731.5 7075799

> `Minus_FRiP`	0.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_gene_sort.bed` (74579,74580)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 74579;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 74580;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_gene_coverage.bed` (74583)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.557GB.  
  PID: 74583;	Command: bedtools;	Return code: 0;	Memory used: 0.061GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed.gz` (74615)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74615;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed` (74616)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 74616;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 12:27:07) elapsed: 70.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed` (74625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 74625;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR"` (74627)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74627;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR_sort.bed` (74628,74629,74630,74631)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 74628;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 74629;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 74631;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB  
  PID: 74630;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_plus_coverage.bed` (74633)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 74633;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_minus_coverage.bed` (74646)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 74646;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR"` (74658)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74658;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR_sort.bed` (74659,74660,74661,74662)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 74659;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 74660;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 74662;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 74661;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_plus_coverage.bed` (74664)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 74664;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_minus_coverage.bed` (74680)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 74680;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer_sort.bed` (74692,74693,74694,74695)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 74692;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 74693;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 74695;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 74694;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_plus_coverage.bed` (74698)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 74698;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_minus_coverage.bed` (74711)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 74711;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon_sort.bed` (74723,74724,74725,74726)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.557GB.  
  PID: 74723;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 74724;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 74726;	Command: bedtools;	Return code: 0;	Memory used: 0.173GB  
  PID: 74725;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_plus_coverage.bed` (74731)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.557GB.  
  PID: 74731;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_minus_coverage.bed` (74747)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.557GB.  
  PID: 74747;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron_sort.bed` (74761,74762,74763,74764)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 74761;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 74763;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 74762;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 74764;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_plus_coverage.bed` (74767)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.557GB.  
  PID: 74767;	Command: bedtools;	Return code: 0;	Memory used: 0.081GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_minus_coverage.bed` (74781)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.557GB.  
  PID: 74781;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_sort.bed` (74796,74797,74798,74799)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74796;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 74797;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 74799;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 74798;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_plus_coverage.bed` (74802)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.557GB.  
  PID: 74802;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_minus_coverage.bed` (74818)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.557GB.  
  PID: 74818;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region"` (74831)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 74831;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed` (74832,74833,74834,74835)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 74832;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 74834;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 74833;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 74835;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (74844)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 74844;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed` (75093)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 75093;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


### Plot cFRiF/FRiF (02-27 12:30:27) elapsed: 201.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_2 -z 3099922541 -n 10264898 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (75117)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.557GB.  
  PID: 75117;	Command: Rscript;	Return code: 0;	Memory used: 0.469GB

> `cFRiF`	QC_hg38/H9_treated_PRO-seq_2_cFRiF.pdf	cFRiF	QC_hg38/H9_treated_PRO-seq_2_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_2 -z 3099922541 -n 10264898 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (75158)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.557GB.  
  PID: 75158;	Command: Rscript;	Return code: 0;	Memory used: 0.466GB

> `FRiF`	QC_hg38/H9_treated_PRO-seq_2_FRiF.pdf	FRiF	QC_hg38/H9_treated_PRO-seq_2_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 12:31:25) elapsed: 57.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_exons_sort.bed` (75189,75190)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.557GB.  
  PID: 75190;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB  
  PID: 75189;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_introns_sort.bed` (75196,75197,75198)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.557GB.  
  PID: 75196;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 75198;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 75197;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_coverage.bed` (75205)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.557GB.  
  PID: 75205;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_coverage.bed` (75231)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.557GB.  
  PID: 75231;	Command: bedtools;	Return code: 0;	Memory used: 0.055GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.8807315)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_rpkm.bed` (75256,75257,75258)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 75256;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 75258;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 75257;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.8807315)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_rpkm.bed` (75261,75262,75263)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 75261;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 75263;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 75262;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed` (75265,75266,75267)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 75265;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 75267;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 75266;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.16	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed --annotate` (75274)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 6.557GB.  
  PID: 75274;	Command: Rscript;	Return code: 0;	Memory used: 0.279GB

> `mRNA contamination`	QC_hg38/H9_treated_PRO-seq_2_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_treated_PRO-seq_2_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed` (75292)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 75292;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (02-27 12:32:36) elapsed: 71.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam` (75300)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.557GB.  
  PID: 75300;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (75310)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_2_plus_cuttrace_vmp2vfp1'
Processing with 2 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270720v1_random', 'chr14_KI270723v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270372v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_exact_body_0-mer.bw'
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:50. Running peak memory: 6.557GB.  
  PID: 75310;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.602GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam` (77100)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.557GB.  
  PID: 77100;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (77107)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_2_minus_cuttrace_9v8bzhmn'
Processing with 2 cores...
Discarding 103 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr22_KI270736v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270757v1']
Keeping 92 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 92 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_exact_body_0-mer.bw'
Merging 92 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:35. Running peak memory: 6.557GB.  
  PID: 77107;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.148GB

Starting cleanup: 77 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:27:15
*  Total elapsed time (all runs):  6:42:17
*         Peak memory (this run):  6.5575 GB
*        Pipeline completed time: 2020-02-27 12:53:19
