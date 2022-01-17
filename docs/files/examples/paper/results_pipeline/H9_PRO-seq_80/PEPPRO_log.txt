### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_80 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_80pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_80pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba26-35c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/
*  Pipeline started at:   (06-15 07:17:27) elapsed: 2.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_80pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_80pct_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_PRO-seq_80`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_80pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_80pct_PE2.fastq.gz

> `File_mb`	1912.05	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:28) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/H9_PRO-seq_80_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_80pct_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/H9_PRO-seq_80_R1.fastq.gz` (217613)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 217613;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/H9_PRO-seq_80_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/H9_PRO-seq_80_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_80pct_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/H9_PRO-seq_80_R2.fastq.gz` (217615)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 217615;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/H9_PRO-seq_80_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/H9_PRO-seq_80_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1.fastq` (217618)
<pre>
</pre>
Command completed. Elapsed time: 0:02:50. Running peak memory: 0.002GB.  
  PID: 217618;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/H9_PRO-seq_80_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2.fastq` (350866)
<pre>
</pre>
Command completed. Elapsed time: 0:02:20. Running peak memory: 0.002GB.  
  PID: 350866;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	92574794	PEPPRO	_RES_

> `Fastq_reads`	92574794	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/H9_PRO-seq_80_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/H9_PRO-seq_80_R2.fastq.gz']

### FASTQ processing:  (06-15 07:25:03) elapsed: 455.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80_R1_cutadapt.txt` (171390)
<pre>
</pre>
Command completed. Elapsed time: 0:02:03. Running peak memory: 3.293GB.  
  PID: 171390;	Command: cutadapt;	Return code: 0;	Memory used: 3.293GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq` (332651,332652)
<pre>
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 3.293GB.  
  PID: 332652;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB  
  PID: 332651;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB

Evaluating read trimming

> `Trimmed_reads`	22665363	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq` (332727)
<pre>
Started analysis of H9_PRO-seq_80_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_80_R1_processed.fastq
Analysis complete for H9_PRO-seq_80_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 3.293GB.  
  PID: 332727;	Command: fastqc;	Return code: 0;	Memory used: 0.181GB

> `FastQC report r1`	fastqc/H9_PRO-seq_80_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_noadap.fastq` (332789)
<pre>
[INFO][0m 2112480 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 3.293GB.  
  PID: 332789;	Command: seqkit;	Return code: 0;	Memory used: 2.037GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_trimmed.fastq` (332903,332931)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.293GB.  
  PID: 332903;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 332931;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	36677941.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	22889260.0	PEPPRO	_RES_

> `Duplicate_reads`	2112480.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.4503	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/processed_R1.flag` (333326)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.293GB.  
  PID: 333326;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80_R2_cutadapt.txt` (333328)
<pre>
</pre>
Command completed. Elapsed time: 0:02:02. Running peak memory: 3.417GB.  
  PID: 333328;	Command: cutadapt;	Return code: 0;	Memory used: 3.417GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed.fastq` (333521,333522)
<pre>
</pre>
Command completed. Elapsed time: 0:01:14. Running peak memory: 3.417GB.  
  PID: 333521;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 333522;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	45330726	PEPPRO	_RES_

> `Trim_loss_rate`	51.03	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq` (333853)
<pre>
Started analysis of H9_PRO-seq_80_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_80_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_80_R1_processed.fastq
Analysis complete for H9_PRO-seq_80_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:49. Running peak memory: 3.417GB.  
  PID: 333853;	Command: fastqc;	Return code: 0;	Memory used: 0.18GB

> `FastQC report r1`	fastqc/H9_PRO-seq_80_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed.fastq` (333926)
<pre>
Started analysis of H9_PRO-seq_80_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_80_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_80_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_80_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:01:00. Running peak memory: 3.417GB.  
  PID: 333926;	Command: fastqc;	Return code: 0;	Memory used: 0.172GB

> `FastQC report r2`	fastqc/H9_PRO-seq_80_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80.histogram`  

> `fastq_pair -t 83317314 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_noadap.fastq` (334042)
<pre>
Left paired: 23236163		Right paired: 23236163
Left single: 161974		Right single: 2141996
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:42. Running peak memory: 4.556GB.  
  PID: 334042;	Command: fastq_pair;	Return code: 0;	Memory used: 4.556GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_noadap.fastq.paired.fq -o H9_PRO-seq_80 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt` (334510)
<pre>
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 4.556GB.  
  PID: 334510;	Command: flash;	Return code: 0;	Memory used: 0.12GB


### Plot adapter insertion distribution (06-15 07:41:58) elapsed: 1015.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt -u 8` (334714)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.556GB.  
  PID: 334714;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_80_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_80_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:42:03) elapsed: 5.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/cutadapt/H9_PRO-seq_80.hist`

> `Degradation_ratio`	1.0312	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed_dups.fastq` (334747)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 4.556GB.  
  PID: 334747;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/processed_R2.flag` (334771)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.556GB.  
  PID: 334771;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/repaired.flag`  

> `fastq_pair -t 83317314 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed.fastq` (334772)
<pre>
Left paired: 22486354		Right paired: 22486354
Left single: 179009		Right single: 2199465
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:18. Running peak memory: 4.823GB.  
  PID: 334772;	Command: fastq_pair;	Return code: 0;	Memory used: 4.823GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq` (335063)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.823GB.  
  PID: 335063;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed.fastq` (335065)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.823GB.  
  PID: 335065;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/repaired.flag` (335066)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.823GB.  
  PID: 335066;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/dups_repaired.flag`  

> `fastq_pair -t 83317314 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed_dups.fastq` (335067)
<pre>
Left paired: 20537390		Right paired: 20537390
Left single: 132325		Right single: 4148429
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:01:49. Running peak memory: 4.894GB.  
  PID: 335067;	Command: fastq_pair;	Return code: 0;	Memory used: 4.894GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_trimmed.fastq` (335354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.894GB.  
  PID: 335354;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed_dups.fastq` (335355)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 335355;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/dups_repaired.flag` (335357)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 335357;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:46:37) elapsed: 274.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:46:37) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/human_rDNA_bt2` (335358)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 335358;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_R2.fq` (335359)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_80 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
22486354 reads; of these:
  22486354 (100.00%) were unpaired; of these:
    20155412 (89.63%) aligned 0 times
    2330942 (10.37%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.37% overall alignment rate

> `Aligned_reads_human_rDNA`	4661884.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.28	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:49:26) elapsed: 168.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/human_rDNA_dups_bt2` (335693)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 335693;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_dups_R2.fq` (335695)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_80 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/fastq/H9_PRO-seq_80_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2330942 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:52:13) elapsed: 167.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_80 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/tmp2gz8zayr -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp.bam` (336112,336113,336114)
<pre>
1845871 reads skipped
0 reads lost
20155412 reads; of these:
  20155412 (100.00%) were paired; of these:
    8354622 (41.45%) aligned concordantly 0 times
    9678102 (48.02%) aligned concordantly exactly 1 time
    2122688 (10.53%) aligned concordantly >1 times
    ----
    8354622 pairs aligned concordantly 0 times; of these:
      1859562 (22.26%) aligned discordantly 1 time
    ----
    6495060 pairs aligned 0 times concordantly or discordantly; of these:
      12990120 mates make up the pairs; of these:
        5331336 (41.04%) aligned 0 times
        2799697 (21.55%) aligned exactly 1 time
        4859087 (37.41%) aligned >1 times
86.77% overall alignment rate
[bam_sort_core] merging from 11 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:39:41. Running peak memory: 4.894GB.  
  PID: 336112;	Command: bowtie2;	Return code: 0;	Memory used: 3.773GB  
  PID: 336113;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 336114;	Command: samtools;	Return code: 0;	Memory used: 0.9GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort.bam` (340181)
<pre>
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 4.894GB.  
  PID: 340181;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	34979488	PEPPRO	_RES_

> `QC_filtered_reads`	20444644	PEPPRO	_RES_

> `Aligned_reads`	14534843.5	PEPPRO	_RES_

> `Alignment_rate`	32.06	PEPPRO	_RES_

> `Total_efficiency`	15.7	PEPPRO	_RES_

> `Read_depth`	2.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_PRO-seq_80 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/tmp2gz8zayr -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp_dups.bam` (341287,341293,341294)
<pre>
18691519 reads; of these:
  18691519 (100.00%) were paired; of these:
    7582143 (40.56%) aligned concordantly 0 times
    9134278 (48.87%) aligned concordantly exactly 1 time
    1975098 (10.57%) aligned concordantly >1 times
    ----
    7582143 pairs aligned concordantly 0 times; of these:
      1760249 (23.22%) aligned discordantly 1 time
    ----
    5821894 pairs aligned 0 times concordantly or discordantly; of these:
      11643788 mates make up the pairs; of these:
        4858981 (41.73%) aligned 0 times
        2643468 (22.70%) aligned exactly 1 time
        4141339 (35.57%) aligned >1 times
87.00% overall alignment rate
[bam_sort_core] merging from 10 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:33:00. Running peak memory: 4.894GB.  
  PID: 341293;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 341287;	Command: bowtie2;	Return code: 0;	Memory used: 3.772GB  
  PID: 341294;	Command: samtools;	Return code: 0;	Memory used: 0.898GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort_dups.bam` (344597)
<pre>
</pre>
Command completed. Elapsed time: 0:01:09. Running peak memory: 4.894GB.  
  PID: 344597;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 09:18:23) elapsed: 5171.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_R2.fq` (344674)
<pre>
</pre>
Command completed. Elapsed time: 0:00:18. Running peak memory: 4.894GB.  
  PID: 344674;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/prealignments/H9_PRO-seq_80_human_rDNA_unmap_R1.fq` (344704)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 4.894GB.  
  PID: 344704;	Command: pigz;	Return code: 0;	Memory used: 0.013GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp.bam` (344738)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 4.894GB.  
  PID: 344738;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	672784	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort.bam` (344767)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 4.894GB.  
  PID: 344767;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/chr_sizes.bed` (344787,344788,344789,344790)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 344788;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 344790;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 344787;	Command: samtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 344789;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_noMT.bam` (344792)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 4.894GB.  
  PID: 344792;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort.bam` (345066)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 345066;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort.bam` (345068)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 4.894GB.  
  PID: 345068;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


### Split BAM file (06-15 09:20:46) elapsed: 143.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam` (345087,345088)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:05. Running peak memory: 4.894GB.  
  PID: 345087;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 345088;	Command: samtools;	Return code: 0;	Memory used: 4.156GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE2.bam` (345232,345233)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:43. Running peak memory: 4.894GB.  
  PID: 345232;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 345233;	Command: samtools;	Return code: 0;	Memory used: 3.33GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp_dups.bam` (345595)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 4.894GB.  
  PID: 345595;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_dups_PE1.bam` (345623,345624)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:00. Running peak memory: 4.894GB.  
  PID: 345623;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 345624;	Command: samtools;	Return code: 0;	Memory used: 3.997GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_dups_PE2.bam` (345768,345769)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:40. Running peak memory: 4.894GB.  
  PID: 345768;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 345769;	Command: samtools;	Return code: 0;	Memory used: 3.21GB


### Calculate library complexity (06-15 09:29:28) elapsed: 522.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_dups_PE1.bam` (345890)
<pre>
BAM_INPUT
TOTAL READS     = 14983475
COUNTS_SUM      = 14983475
DISTINCT READS  = 1.23124e+07
DISTINCT COUNTS = 284
MAX COUNT       = 22754
COUNTS OF 1     = 1.10363e+07
OBSERVED COUNTS (22755)
1	11036268
2	873369
3	200465
4	78869
5	39923
6	23143
7	14523
8	9658
9	6943
10	5057
11	3773
12	2955
13	2220
14	1753
15	1498
16	1200
17	1075
18	884
19	794
20	658
21	562
22	518
23	441
24	404
25	329
26	321
27	289
28	282
29	236
30	224
31	176
32	185
33	177
34	132
35	123
36	153
37	135
38	109
39	113
40	109
41	95
42	88
43	94
44	75
45	90
46	67
47	79
48	74
49	61
50	66
51	52
52	55
53	44
54	39
55	43
56	49
57	34
58	37
59	46
60	34
61	29
62	24
63	37
64	22
65	26
66	33
67	27
68	34
69	26
70	26
71	26
72	19
73	22
74	16
75	18
76	21
77	20
78	16
79	21
80	13
81	22
82	17
83	14
84	14
85	13
86	13
87	5
88	11
89	14
90	10
91	6
92	11
93	13
94	14
95	6
96	9
97	9
98	10
99	12
100	7
101	11
102	7
103	8
104	3
105	7
106	4
107	5
108	10
109	11
110	8
111	6
112	10
113	4
114	8
115	6
116	7
117	4
118	4
119	5
120	2
121	10
122	2
123	5
124	3
125	3
126	3
127	7
128	7
129	6
130	4
131	5
132	3
133	2
134	2
135	2
136	2
137	6
138	6
139	2
140	2
141	5
142	2
143	2
144	3
145	6
146	4
147	3
148	3
149	1
153	1
154	4
155	1
156	2
157	1
158	2
159	1
160	2
161	4
162	2
163	2
164	7
165	7
166	3
168	1
169	3
170	3
171	1
172	3
173	1
174	1
175	2
177	2
178	1
179	3
180	2
181	1
182	2
183	2
184	4
186	1
188	6
189	3
192	1
194	1
195	2
196	1
198	2
201	1
203	2
205	1
206	1
207	1
209	1
210	1
214	2
215	1
217	2
218	1
219	2
220	1
222	2
223	1
224	1
226	1
227	1
228	1
230	1
232	1
233	1
243	1
244	1
248	2
252	2
253	1
254	2
256	1
257	2
259	3
261	1
264	1
265	2
269	2
272	1
274	1
275	2
277	1
281	1
282	1
285	1
292	1
294	2
303	1
304	1
307	1
309	1
311	1
315	1
316	1
327	1
334	1
346	1
351	1
355	1
361	1
378	1
395	1
397	1
406	1
410	1
412	1
418	1
420	1
430	1
443	1
446	1
479	1
555	1
567	1
569	2
571	1
579	1
583	1
716	1
719	1
741	1
786	1
815	1
956	1
1013	1
1067	1
1304	1
1532	1
1541	1
1558	1
1869	1
1983	1
2093	1
2831	1
3629	1
7419	1
8643	1
11051	1
20426	1
22754	1

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
Command completed. Elapsed time: 0:01:23. Running peak memory: 4.894GB.  
  PID: 345890;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_dups_PE1.bam` (346215)
<pre>
BAM_INPUT
TOTAL READS     = 14983475
DISTINCT READS  = 1.23124e+07
DISTINCT COUNTS = 284
MAX COUNT       = 22754
COUNTS OF 1     = 1.10363e+07
MAX TERMS       = 100
OBSERVED COUNTS (22755)
1	11036268
2	873369
3	200465
4	78869
5	39923
6	23143
7	14523
8	9658
9	6943
10	5057
11	3773
12	2955
13	2220
14	1753
15	1498
16	1200
17	1075
18	884
19	794
20	658
21	562
22	518
23	441
24	404
25	329
26	321
27	289
28	282
29	236
30	224
31	176
32	185
33	177
34	132
35	123
36	153
37	135
38	109
39	113
40	109
41	95
42	88
43	94
44	75
45	90
46	67
47	79
48	74
49	61
50	66
51	52
52	55
53	44
54	39
55	43
56	49
57	34
58	37
59	46
60	34
61	29
62	24
63	37
64	22
65	26
66	33
67	27
68	34
69	26
70	26
71	26
72	19
73	22
74	16
75	18
76	21
77	20
78	16
79	21
80	13
81	22
82	17
83	14
84	14
85	13
86	13
87	5
88	11
89	14
90	10
91	6
92	11
93	13
94	14
95	6
96	9
97	9
98	10
99	12
100	7
101	11
102	7
103	8
104	3
105	7
106	4
107	5
108	10
109	11
110	8
111	6
112	10
113	4
114	8
115	6
116	7
117	4
118	4
119	5
120	2
121	10
122	2
123	5
124	3
125	3
126	3
127	7
128	7
129	6
130	4
131	5
132	3
133	2
134	2
135	2
136	2
137	6
138	6
139	2
140	2
141	5
142	2
143	2
144	3
145	6
146	4
147	3
148	3
149	1
153	1
154	4
155	1
156	2
157	1
158	2
159	1
160	2
161	4
162	2
163	2
164	7
165	7
166	3
168	1
169	3
170	3
171	1
172	3
173	1
174	1
175	2
177	2
178	1
179	3
180	2
181	1
182	2
183	2
184	4
186	1
188	6
189	3
192	1
194	1
195	2
196	1
198	2
201	1
203	2
205	1
206	1
207	1
209	1
210	1
214	2
215	1
217	2
218	1
219	2
220	1
222	2
223	1
224	1
226	1
227	1
228	1
230	1
232	1
233	1
243	1
244	1
248	2
252	2
253	1
254	2
256	1
257	2
259	3
261	1
264	1
265	2
269	2
272	1
274	1
275	2
277	1
281	1
282	1
285	1
292	1
294	2
303	1
304	1
307	1
309	1
311	1
315	1
316	1
327	1
334	1
346	1
351	1
355	1
361	1
378	1
395	1
397	1
406	1
410	1
412	1
418	1
420	1
430	1
443	1
446	1
479	1
555	1
567	1
569	2
571	1
579	1
583	1
716	1
719	1
741	1
786	1
815	1
956	1
1013	1
1067	1
1304	1
1532	1
1541	1
1558	1
1869	1
1983	1
2093	1
2831	1
3629	1
7419	1
8643	1
11051	1
20426	1
22754	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
..........................._...................................._.....................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:30. Running peak memory: 4.894GB.  
  PID: 346215;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_counts.txt` (346289)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 4.894GB.  
  PID: 346289;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_plot` (346314)
<pre>
Processing H9_PRO-seq_80
INFO: Found real counts for H9_PRO-seq_80 - Total (M): 15.588108 Unique (M): 14.983475

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.894GB.  
  PID: 346314;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Library complexity`	QC_hg38/H9_PRO-seq_80_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_80_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8537	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 09:32:48) elapsed: 200.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam` (346333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 4.894GB.  
  PID: 346333;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_bamQC.tsv` (346344)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/tmp_H9_PRO-seq_80_PE1_5q6gvnk9'
Processing with 12 cores...
Discarding 101 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270755v1']
Keeping 94 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 4.894GB.  
  PID: 346344;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.548GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	7794054.0	PEPPRO	_RES_

> `PBC2`	7794054.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_unmap.bam` (346385)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 4.894GB.  
  PID: 346385;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_temp.bam`

> `Unmapped_reads`	5331336	PEPPRO	_RES_

### Split BAM by strand (06-15 09:33:25) elapsed: 37.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam` (346428)
<pre>
</pre>
Command completed. Elapsed time: 0:00:54. Running peak memory: 4.894GB.  
  PID: 346428;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam` (346474)
<pre>
</pre>
Command completed. Elapsed time: 0:00:52. Running peak memory: 4.894GB.  
  PID: 346474;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 09:35:11) elapsed: 106.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (346707)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 346707;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_plus_TssEnrichment.txt` (346708)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.894GB.  
  PID: 346708;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.974GB


> `TSS_coding_score`	33.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_minus_TssEnrichment.txt` (346739)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 4.894GB.  
  PID: 346739;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.992GB


> `TSS_non-coding_score`	12.3	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_minus_TssEnrichment.txt` (346770)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.894GB.  
  PID: 346770;	Command: Rscript;	Return code: 0;	Memory used: 0.321GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_80_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_80_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt` (346792,346793,346794,346795)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 346792;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 346794;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 346793;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 346795;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_keep.txt` (346797)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 346797;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 09:35:27) elapsed: 16.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_ensembl_tss.bed` (346799,346800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 4.894GB.  
  PID: 346799;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 346800;	Command: bedtools;	Return code: 0;	Memory used: 0.098GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_ensembl_gene_body.bed` (346804,346805)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 346804;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 346805;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_TSS_density.bed` (346807,346808,346809,346810)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 4.894GB.  
  PID: 346807;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 346809;	Command: sort;	Return code: 0;	Memory used: 0.012GB  
  PID: 346808;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 346810;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_gene_body_density.bed` (346830,346831,346832)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 4.894GB.  
  PID: 346831;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 346830;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB  
  PID: 346832;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/tmp3j3go7g0` (346854,346855,346856)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 346854;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 346856;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 346855;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/tmp3j3go7g0 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0155008) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/tmp3j3go7g0 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_pause_index.bed` (346862)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 346862;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	33.72	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_pause_index.bed` (346867)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.894GB.  
  PID: 346867;	Command: Rscript;	Return code: 0;	Memory used: 0.319GB

> `Pause index`	QC_hg38/H9_PRO-seq_80_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_80_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_pause_index.bed` (346888)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 346888;	Command: pigz;	Return code: 0;	Memory used: 0.004GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 09:36:19) elapsed: 52.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam`
14534843.5 5639799

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam`
14534843.5 5259266

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_gene_sort.bed` (346927,346928)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.894GB.  
  PID: 346927;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 346928;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_gene_coverage.bed` (346931)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 4.894GB.  
  PID: 346931;	Command: bedtools;	Return code: 0;	Memory used: 0.085GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/hg38_annotations.bed.gz` (346969)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 346969;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/hg38_annotations.bed` (346970)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 346970;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 09:37:11) elapsed: 53.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/raw/hg38_annotations.bed` (346978)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.894GB.  
  PID: 346978;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Enhancer_sort.bed` (346981,346982,346983,346984)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.894GB.  
  PID: 346981;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 346982;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 346984;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB  
  PID: 346983;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Enhancer_plus_coverage.bed` (346987)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 4.894GB.  
  PID: 346987;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Enhancer_minus_coverage.bed` (346997)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 4.894GB.  
  PID: 346997;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_sort.bed` (347007,347008,347009,347010)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 347007;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 347008;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 347010;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 347009;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_plus_coverage.bed` (347013)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.894GB.  
  PID: 347013;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_minus_coverage.bed` (347024)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 4.894GB.  
  PID: 347024;	Command: bedtools;	Return code: 0;	Memory used: 0.016GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_Flanking_Region"` (347034)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 347034;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_Flanking_Region_sort.bed` (347035,347036,347037,347038)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.894GB.  
  PID: 347035;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 347037;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 347036;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 347038;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_Flanking_Region_plus_coverage.bed` (347041)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.894GB.  
  PID: 347041;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_Flanking_Region_minus_coverage.bed` (347051)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 4.894GB.  
  PID: 347051;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/5_UTR"` (347061)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 347061;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/5_UTR_sort.bed` (347062,347063,347064,347065)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.894GB.  
  PID: 347062;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 347063;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 347065;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 347064;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_5_UTR_plus_coverage.bed` (347068)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.894GB.  
  PID: 347068;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_5_UTR_minus_coverage.bed` (347095)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.894GB.  
  PID: 347095;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/3_UTR"` (347107)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 347107;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/3_UTR_sort.bed` (347108,347109,347110,347111)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.894GB.  
  PID: 347108;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 347109;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 347111;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 347110;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_3_UTR_plus_coverage.bed` (347128)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.894GB.  
  PID: 347128;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_3_UTR_minus_coverage.bed` (347140)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 4.894GB.  
  PID: 347140;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Exon_sort.bed` (347150,347151,347152,347153)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 4.894GB.  
  PID: 347150;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 347151;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 347153;	Command: bedtools;	Return code: 0;	Memory used: 0.171GB  
  PID: 347152;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Exon_plus_coverage.bed` (347159)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.894GB.  
  PID: 347159;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Exon_minus_coverage.bed` (347204)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.894GB.  
  PID: 347204;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Intron_sort.bed` (347216,347217,347218,347219)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.894GB.  
  PID: 347216;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 347218;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 347217;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 347219;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Intron_plus_coverage.bed` (347222)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 4.894GB.  
  PID: 347222;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Intron_minus_coverage.bed` (347239)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 4.894GB.  
  PID: 347239;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB


### Plot cFRiF/FRiF (06-15 09:39:52) elapsed: 161.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_80 -z 3099922541 -n 7955869 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Intron_plus_coverage.bed` (347266)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 4.894GB.  
  PID: 347266;	Command: Rscript;	Return code: 0;	Memory used: 0.484GB

> `cFRiF`	QC_hg38/H9_PRO-seq_80_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_80_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_80 -z 3099922541 -n 7955869 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_Intron_plus_coverage.bed` (347560)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 4.894GB.  
  PID: 347560;	Command: Rscript;	Return code: 0;	Memory used: 0.484GB

> `FRiF`	QC_hg38/H9_PRO-seq_80_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_80_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 09:40:52) elapsed: 60.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_exons_sort.bed` (347596,347597)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.894GB.  
  PID: 347597;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB  
  PID: 347596;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_introns_sort.bed` (347619,347620,347621)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.894GB.  
  PID: 347619;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 347621;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 347620;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exons_coverage.bed` (347627)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 4.894GB.  
  PID: 347627;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_introns_coverage.bed` (347645)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 4.894GB.  
  PID: 347645;	Command: bedtools;	Return code: 0;	Memory used: 0.073GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.5348435)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exons_rpkm.bed` (347665,347666,347667)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.894GB.  
  PID: 347665;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 347667;	Command: sort;	Return code: 0;	Memory used: 0.007GB  
  PID: 347666;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/14.5348435)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_introns_rpkm.bed` (347669,347670,347671)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.894GB.  
  PID: 347669;	Command: awk;	Return code: 0;	Memory used: 0.007GB  
  PID: 347671;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 347670;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exon_intron_ratios.bed` (347674,347675,347676)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 347674;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 347676;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 347675;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.29	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exon_intron_ratios.bed --annotate` (347682)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 4.894GB.  
  PID: 347682;	Command: Rscript;	Return code: 0;	Memory used: 0.321GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_80_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_80_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/QC_hg38/H9_PRO-seq_80_exon_intron_ratios.bed` (347702)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.894GB.  
  PID: 347702;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (06-15 09:41:50) elapsed: 58.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam` (347711)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.894GB.  
  PID: 347711;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 14534843.5` (347716)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_80_plus_cuttrace_konb2rr9'
Processing with 4 cores...
stdin is empty of data
Discarding 113 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_GL000009v2_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 82 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_plus_exact_body_0-mer.bw'
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:16. Running peak memory: 4.894GB.  
  PID: 347716;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.768GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam` (349251)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 4.894GB.  
  PID: 349251;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 14534843.5` (349257)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/aligned_hg38/H9_PRO-seq_80_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_80_minus_cuttrace_tx0_9dvj'
Processing with 4 cores...
stdin is empty of data
Discarding 113 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr14_GL000009v2_random', 'chr14_KI270723v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270755v1']
Keeping 82 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270731v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270509v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270362v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_minus_exact_body_0-mer.bw'
Merging 82 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_80/signal_hg38/H9_PRO-seq_80_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:11. Running peak memory: 4.894GB.  
  PID: 349257;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.34GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:39:05
*  Total elapsed time (all runs):  4:36:50
*         Peak memory (this run):  4.8945 GB
*        Pipeline completed time: 2020-06-15 09:56:30
