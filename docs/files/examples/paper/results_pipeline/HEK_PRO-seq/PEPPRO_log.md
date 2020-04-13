### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name HEK_PRO-seq --genome hg38 --input /project/shefflab/data/guertin/fastq/SRR8608074_PE1.fastq.gz --single-or-paired paired --input2 /project/shefflab/data/guertin/fastq/SRR8608074_PE2.fastq.gz --protocol PRO --umi-len 8 --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 8 -M 16000`
*         Compute host:  udc-aj38-13c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/
*  Pipeline started at:   (02-27 09:25:04) elapsed: 2.0 _TIME_

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
*              `input`:  `['/project/shefflab/data/guertin/fastq/SRR8608074_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data/guertin/fastq/SRR8608074_PE2.fastq.gz']`
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
*        `sample_name`:  `HEK_PRO-seq`
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

Local input file: /project/shefflab/data/guertin/fastq/SRR8608074_PE1.fastq.gz
Local input file: /project/shefflab/data/guertin/fastq/SRR8608074_PE2.fastq.gz

> `File_mb`	1337.94	PEPPRO	_RES_

> `Read_type`	paired	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-27 09:25:05) elapsed: 0.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R1.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/SRR8608074_PE1.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R1.fastq.gz` (184008)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 184008;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R2.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/SRR8608074_PE2.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R2.fastq.gz` (184009)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 184009;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1.fastq`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2.fastq`  

> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R1.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1.fastq` (184010)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 0.003GB.  
  PID: 184010;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R2.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2.fastq` (184052)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 0.003GB.  
  PID: 184052;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	51533198	PEPPRO	_RES_

> `Fastq_reads`	51533198	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/HEK_PRO-seq_R2.fastq.gz']

### FASTQ processing:  (02-27 09:26:37) elapsed: 93.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R1_cutadapt.txt` (184144)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.276GB.  
  PID: 184144;	Command: cutadapt;	Return code: 0;	Memory used: 0.276GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq` (184204,184205)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 0.276GB.  
  PID: 184204;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 184205;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	20103582	PEPPRO	_RES_

> `Trim_loss_rate`	60.99	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq` (184240)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
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
Command completed. Elapsed time: 0:00:49. Running peak memory: 0.276GB.  
  PID: 184240;	Command: fastqc;	Return code: 0;	Memory used: 0.162GB

> `FastQC report r1`	fastqc/HEK_PRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq`  

> `seqkit rmdup --threads 8 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_dedup.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq` (184316)
<pre>
[INFO][0m 916082 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 0.991GB.  
  PID: 184316;	Command: seqkit;	Return code: 0;	Memory used: 0.991GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq` (184375,184376)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 0.991GB.  
  PID: 184375;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 184376;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	16626152.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	5380555.0	PEPPRO	_RES_

> `Duplicate_reads`	916082.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	20.8819	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/processed_R1.flag` (184422)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.991GB.  
  PID: 184422;	Command: touch;	Return code: 0;	Memory used: 0.002GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq`  

> `(cutadapt -j 8 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_R2_cutadapt.txt` (184424)
<pre>
</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 0.991GB.  
  PID: 184424;	Command: cutadapt;	Return code: 0;	Memory used: 0.283GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq` (184478,184479)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 0.991GB.  
  PID: 184478;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 184479;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB

Evaluating read trimming

> `Trimmed_reads`	40207164	PEPPRO	_RES_

> `Trim_loss_rate`	21.98	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq` (184733)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
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
Command completed. Elapsed time: 0:00:53. Running peak memory: 0.991GB.  
  PID: 184733;	Command: fastqc;	Return code: 0;	Memory used: 0.162GB

> `FastQC report r1`	fastqc/HEK_PRO-seq_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq` (184796)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
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
Command completed. Elapsed time: 0:00:48. Running peak memory: 0.991GB.  
  PID: 184796;	Command: fastqc;	Return code: 0;	Memory used: 0.189GB

> `FastQC report r2`	fastqc/HEK_PRO-seq_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.histogram`  

> `fastq_pair -t 46379878 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq` (184858)
<pre>
Left paired: 20354925		Right paired: 20354925
Left single: 31120		Right single: 456238
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:27. Running peak memory: 2.47GB.  
  PID: 184858;	Command: fastq_pair;	Return code: 0;	Memory used: 2.47GB


> `flash -q -t 8 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_noadap.fastq.paired.fq -o HEK_PRO-seq -d /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt` (184938)
<pre>
</pre>
Command completed. Elapsed time: 0:01:02. Running peak memory: 2.47GB.  
  PID: 184938;	Command: flash;	Return code: 0;	Memory used: 0.12GB


### Plot adapter insertion distribution (02-27 09:34:34) elapsed: 477.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt -u 8` (185146)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 2.47GB.  
  PID: 185146;	Command: Rscript;	Return code: 0;	Memory used: 0.127GB

> `Adapter insertion distribution`	cutadapt/HEK_PRO-seq_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/HEK_PRO-seq_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	37	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:34:42) elapsed: 8.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/cutadapt/HEK_PRO-seq.hist`

> `Degradation_ratio`	0.5783	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq` (185175)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 2.47GB.  
  PID: 185175;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/processed_R2.flag` (185188)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.47GB.  
  PID: 185188;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/repaired.flag`  

> `fastq_pair -t 46379878 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq` (185190)
<pre>
Left paired: 20071296		Right paired: 20071296
Left single: 32287		Right single: 474552
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:22. Running peak memory: 2.502GB.  
  PID: 185190;	Command: fastq_pair;	Return code: 0;	Memory used: 2.502GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq` (185512)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 2.502GB.  
  PID: 185512;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq` (185514)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.502GB.  
  PID: 185514;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/repaired.flag` (185515)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.502GB.  
  PID: 185515;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/dups_repaired.flag`  

> `fastq_pair -t 46379878 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq` (185516)
<pre>
Left paired: 19171561		Right paired: 19171561
Left single: 25508		Right single: 1374287
Writing the paired reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:20. Running peak memory: 2.591GB.  
  PID: 185516;	Command: fastq_pair;	Return code: 0;	Memory used: 2.591GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq` (185617)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 2.591GB.  
  PID: 185617;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq` (185618)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 2.591GB.  
  PID: 185618;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/dups_repaired.flag` (185620)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.591GB.  
  PID: 185620;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (02-27 09:37:42) elapsed: 180.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:37:42) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_bt2` (185621)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.591GB.  
  PID: 185621;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq` (185622)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HEK_PRO-seq -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
Missing stat 'Aligned_reads_human_rDNA'
20071296 reads; of these:
  20071296 (100.00%) were unpaired; of these:
    17956185 (89.46%) aligned 0 times
    2115111 (10.54%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.54% overall alignment rate

> `Aligned_reads_human_rDNA`	4230222.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.52	PEPPRO	_RES_

### Map to human_rDNA (02-27 09:39:55) elapsed: 133.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_dups_bt2` (185754)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 2.591GB.  
  PID: 185754;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio/peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_dups_R2.fq` (185755)
<pre>
</pre>

> `(bowtie2 -p 8 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id HEK_PRO-seq -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/fastq/HEK_PRO-seq_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2115112 reads skipped
0 reads lost

### Map to genome (02-27 09:42:24) elapsed: 149.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id HEK_PRO-seq -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/tmpviac0ifo -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam` (186114,186115,186116)
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
[bam_sort_core] merging from 8 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:49:27. Running peak memory: 3.659GB.  
  PID: 186115;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 186114;	Command: bowtie2;	Return code: 0;	Memory used: 3.659GB  
  PID: 186116;	Command: samtools;	Return code: 0;	Memory used: 0.884GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam` (191321)
<pre>
</pre>
Command completed. Elapsed time: 0:00:53. Running peak memory: 3.659GB.  
  PID: 191321;	Command: samtools;	Return code: 0;	Memory used: 0.015GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	33104378	PEPPRO	_RES_

> `QC_filtered_reads`	18769818	PEPPRO	_RES_

> `Aligned_reads`	14334559.5	PEPPRO	_RES_

> `Alignment_rate`	35.65	PEPPRO	_RES_

> `Total_efficiency`	27.82	PEPPRO	_RES_

> `Read_depth`	2.41	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam`  

> `bowtie2 -p 8 --very-sensitive -X 2000 --rg-id HEK_PRO-seq -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/tmpviac0ifo -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp_dups.bam` (192829,192830,192831)
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
[bam_sort_core] merging from 8 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:44:59. Running peak memory: 3.663GB.  
  PID: 192829;	Command: bowtie2;	Return code: 0;	Memory used: 3.663GB  
  PID: 192830;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 192831;	Command: samtools;	Return code: 0;	Memory used: 0.882GB


> `samtools view -q 10 -b -@ 8 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp_dups.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam` (197499)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 3.663GB.  
  PID: 197499;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


### Compress all unmapped read files (02-27 11:31:40) elapsed: 6556.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R1.fq` (197567)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.663GB.  
  PID: 197567;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/prealignments/HEK_PRO-seq_human_rDNA_unmap_R2.fq` (197600)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 3.663GB.  
  PID: 197600;	Command: pigz;	Return code: 0;	Memory used: 0.008GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam` (197639)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.663GB.  
  PID: 197639;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	325556	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam` (197682)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.663GB.  
  PID: 197682;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_noMT.bam` (197700,197701,197702,197703)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 3.663GB.  
  PID: 197702;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 197700;	Command: samtools;	Return code: 0;	Memory used: 0.005GB  
  PID: 197701;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 197703;	Command: xargs;	Return code: 0;	Memory used: 0.074GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam` (197741)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 197741;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam` (197743)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.663GB.  
  PID: 197743;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


### Split BAM file (02-27 11:34:08) elapsed: 148.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam` (197761,197762)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:54. Running peak memory: 3.663GB.  
  PID: 197761;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 197762;	Command: samtools;	Return code: 0;	Memory used: 3.277GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE2.bam` (198151,198152)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:41. Running peak memory: 3.663GB.  
  PID: 198151;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 198152;	Command: samtools;	Return code: 0;	Memory used: 2.752GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp_dups.bam` (198327)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 3.663GB.  
  PID: 198327;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE1.bam` (198358,198359)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:53. Running peak memory: 3.663GB.  
  PID: 198358;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 198359;	Command: samtools;	Return code: 0;	Memory used: 3.226GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_sort_dups.bam | samtools sort - -@ 8 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE2.bam` (198699,198700)
<pre>
[bam_sort_core] merging from 0 files and 8 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:39. Running peak memory: 3.663GB.  
  PID: 198699;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 198700;	Command: samtools;	Return code: 0;	Memory used: 2.708GB


### Calculate library complexity (02-27 11:42:15) elapsed: 487.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_out.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE1.bam` (198812)
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
Command completed. Elapsed time: 0:01:27. Running peak memory: 3.663GB.  
  PID: 198812;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE1.bam` (198916)
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
................................................................._...................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:32. Running peak memory: 3.663GB.  
  PID: 198916;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_counts.txt` (199217)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.663GB.  
  PID: 199217;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_plot` (199243)
<pre>
Processing HEK_PRO-seq
INFO: Found real counts for HEK_PRO-seq - Total (M): 15.183469 Unique (M): 14.923343

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.663GB.  
  PID: 199243;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `Library complexity`	QC_hg38/HEK_PRO-seq_preseq_plot.pdf	Library complexity	QC_hg38/HEK_PRO-seq_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.9097	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 11:45:43) elapsed: 207.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam` (199263)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.663GB.  
  PID: 199263;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -c 8 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_bamQC.tsv` (199273)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/tmp_HEK_PRO-seq_PE1_begp3fx0'
Processing with 8 cores...
Discarding 105 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270748v1', 'chrUn_KI270752v1', 'chrUn_KI270756v1']
Keeping 90 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.663GB.  
  PID: 199273;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.084GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	7591734.5	PEPPRO	_RES_

> `PBC2`	7591734.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_unmap.bam`  

> `samtools view -b -@ 8 -f 12  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_unmap.bam` (199310)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.663GB.  
  PID: 199310;	Command: samtools;	Return code: 0;	Memory used: 0.008GB


> `samtools view -c -f 4 -@ 8 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_temp.bam`

> `Unmapped_reads`	2807990	PEPPRO	_RES_

### Split BAM by strand (02-27 11:46:24) elapsed: 42.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam` (199344)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 3.663GB.  
  PID: 199344;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam` (199391)
<pre>
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 3.663GB.  
  PID: 199391;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 11:48:05) elapsed: 101.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (199448)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 199448;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/plus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_plus_TssEnrichment.txt` (199449)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.663GB.  
  PID: 199449;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.611GB


> `TSS_coding_score`	16.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/minus_TSS.tsv -p ends -c 8 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_minus_TssEnrichment.txt` (199474)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.663GB.  
  PID: 199474;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.619GB


> `TSS_non-coding_score`	2.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_minus_TssEnrichment.txt` (199498)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.663GB.  
  PID: 199498;	Command: Rscript;	Return code: 0;	Memory used: 0.238GB

> `TSS enrichment`	QC_hg38/HEK_PRO-seq_TSSenrichment.pdf	TSS enrichment	QC_hg38/HEK_PRO-seq_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt` (199517,199518,199519,199520)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 199517;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 199519;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 199518;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 199520;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt` (199522)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 199522;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (02-27 11:48:24) elapsed: 19.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_tss.bed` (199524,199525)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.663GB.  
  PID: 199524;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 199525;	Command: bedtools;	Return code: 0;	Memory used: 0.1GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed` (199529,199530)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 199529;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 199530;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_TSS_density.bed` (199535,199536,199537,199538)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 3.663GB.  
  PID: 199535;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 199537;	Command: sort;	Return code: 0;	Memory used: 0.014GB  
  PID: 199536;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 199538;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_gene_body_density.bed` (199561,199562,199563)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 3.663GB.  
  PID: 199562;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 199561;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB  
  PID: 199563;	Command: sort;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed` (199587,199588,199589)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 199587;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 199589;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 199588;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	14.83	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed` (199595)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.663GB.  
  PID: 199595;	Command: Rscript;	Return code: 0;	Memory used: 0.236GB

> `Pause index`	QC_hg38/HEK_PRO-seq_pause_index.pdf	Pause index	QC_hg38/HEK_PRO-seq_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_pause_index.bed` (199614)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 199614;	Command: pigz;	Return code: 0;	Memory used: 0.006GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 11:49:22) elapsed: 58.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam`
14334559.5 5618447

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam`
14334559.5 5337871

> `Minus_FRiP`	0.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_gene_sort.bed` (199662,199663)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 199662;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 199663;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_gene_coverage.bed` (199666)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.663GB.  
  PID: 199666;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed.gz` (199904)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 199904;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 8 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed` (199905)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 199905;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 11:50:19) elapsed: 57.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/raw/hg38_annotations.bed` (199914)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 199914;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR"` (199916)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 199916;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR_sort.bed` (199917,199918,199919,199920)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 199917;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 199918;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 199920;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 199919;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_plus_coverage.bed` (199923)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.663GB.  
  PID: 199923;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_minus_coverage.bed` (199948)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.663GB.  
  PID: 199948;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR"` (199958)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 199958;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR_sort.bed` (199959,199960,199961,199962)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 199959;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 199960;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 199962;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 199961;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_plus_coverage.bed` (199964)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.663GB.  
  PID: 199964;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_minus_coverage.bed` (199975)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.663GB.  
  PID: 199975;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer_sort.bed` (199984,199985,199986,199987)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 199984;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 199985;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 199987;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB  
  PID: 199986;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_plus_coverage.bed` (199990)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.663GB.  
  PID: 199990;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_minus_coverage.bed` (200012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.663GB.  
  PID: 200012;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Exon_sort.bed` (200022,200023,200024,200025)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.663GB.  
  PID: 200022;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 200023;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 200025;	Command: bedtools;	Return code: 0;	Memory used: 0.178GB  
  PID: 200024;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_plus_coverage.bed` (200033)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.663GB.  
  PID: 200033;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_minus_coverage.bed` (200045)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.663GB.  
  PID: 200045;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Intron_sort.bed` (200056,200057,200058,200059)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 200056;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 200058;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 200057;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 200059;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_plus_coverage.bed` (200062)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.663GB.  
  PID: 200062;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_minus_coverage.bed` (200086)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.663GB.  
  PID: 200086;	Command: bedtools;	Return code: 0;	Memory used: 0.018GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_sort.bed` (200098,200099,200100,200102)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 200098;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 200099;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 200102;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 200100;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_plus_coverage.bed` (200104)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.663GB.  
  PID: 200104;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_minus_coverage.bed` (200117)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.663GB.  
  PID: 200117;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region"` (200129)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 200129;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed` (200130,200131,200132,200133)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 200130;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 200132;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 200131;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 200133;	Command: bedtools;	Return code: 0;	Memory used: 0.049GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (200136)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.663GB.  
  PID: 200136;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_minus_coverage.bed` (200158)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.663GB.  
  PID: 200158;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


### Plot cFRiF/FRiF (02-27 11:53:01) elapsed: 161.0 _TIME_


> `samtools view -@ 8 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s HEK_PRO-seq -z 3099922541 -n 7704166 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (200180)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 3.663GB.  
  PID: 200180;	Command: Rscript;	Return code: 0;	Memory used: 0.469GB

> `cFRiF`	QC_hg38/HEK_PRO-seq_cFRiF.pdf	cFRiF	QC_hg38/HEK_PRO-seq_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s HEK_PRO-seq -z 3099922541 -n 7704166 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_Promoter_Flanking_Region_plus_coverage.bed` (200225)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 3.663GB.  
  PID: 200225;	Command: Rscript;	Return code: 0;	Memory used: 0.465GB

> `FRiF`	QC_hg38/HEK_PRO-seq_FRiF.pdf	FRiF	QC_hg38/HEK_PRO-seq_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 11:54:07) elapsed: 66.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_exons_sort.bed` (200260,200261)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.663GB.  
  PID: 200261;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 200260;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_introns_sort.bed` (200268,200269,200270)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.663GB.  
  PID: 200268;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 200270;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 200269;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_coverage.bed` (200277)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 3.663GB.  
  PID: 200277;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_PE1.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_coverage.bed` (200299)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 3.663GB.  
  PID: 200299;	Command: bedtools;	Return code: 0;	Memory used: 0.026GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.3345595)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_rpkm.bed` (200379,200386,200387)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 200379;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 200387;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 200386;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.3345595)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_rpkm.bed` (200499,200500,200501)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.663GB.  
  PID: 200499;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 200501;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 200500;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed` (200522,200523,200524)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 200522;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 200524;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 200523;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.25	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed --annotate` (200530)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.663GB.  
  PID: 200530;	Command: Rscript;	Return code: 0;	Memory used: 0.237GB

> `mRNA contamination`	QC_hg38/HEK_PRO-seq_mRNA_contamination.pdf	mRNA contamination	QC_hg38/HEK_PRO-seq_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed.gz`  

> `pigz -f -p 8 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/QC_hg38/HEK_PRO-seq_exon_intron_ratios.bed` (200549)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.663GB.  
  PID: 200549;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (02-27 11:55:08) elapsed: 61.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam` (200558)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.663GB.  
  PID: 200558;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (200567)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_plus.bam'
Temporary files will be stored in: 'tmp_HEK_PRO-seq_plus_cuttrace_wn79bgez'
Processing with 2 cores...
Discarding 121 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1']
Keeping 74 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 74 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_exact_body_0-mer.bw'
Merging 74 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:07. Running peak memory: 3.663GB.  
  PID: 200567;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.53GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam` (202282)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.663GB.  
  PID: 202282;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_smooth_body_0-mer.bw -p 5 --variable-step --tail-edge` (202291)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/aligned_hg38/HEK_PRO-seq_minus.bam'
Temporary files will be stored in: 'tmp_HEK_PRO-seq_minus_cuttrace_jtw4ag67'
Processing with 2 cores...
Discarding 117 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr5_GL000208v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr17_KI270730v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 78 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 78 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_exact_body_0-mer.bw'
Merging 78 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/HEK_PRO-seq/signal_hg38/HEK_PRO-seq_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:03. Running peak memory: 3.663GB.  
  PID: 202291;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.539GB

Starting cleanup: 77 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:50:30
*  Total elapsed time (all runs):  5:35:59
*         Peak memory (this run):  3.6634 GB
*        Pipeline completed time: 2020-02-27 12:15:32
