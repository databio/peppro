### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_20 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 4 -M 8000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-aw29-25a
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/
*  Pipeline started at:   (06-11 17:54:22) elapsed: 5.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.8
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16373 insertions(+), 3522 deletions(-)

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `cutadapt`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `4`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE2.fastq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `8000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `H9_PRO-seq_20`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE2.fastq.gz

> `File_mb`	492.8	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:54:24) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz` (51553)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 51553;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz` (51555)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 51555;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1.fastq` (51556)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 0.002GB.  
  PID: 51556;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2.fastq` (52032)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 0.002GB.  
  PID: 52032;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	23139336	PEPPRO	_RES_

> `Fastq_reads`	23139336	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz']

### FASTQ processing:  (06-11 17:55:14) elapsed: 51.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq`  

> `(cutadapt -j 4 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20_R1_cutadapt.txt` (52485)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 0.923GB.  
  PID: 52485;	Command: cutadapt;	Return code: 0;	Memory used: 0.923GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq` (52688,52690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.923GB.  
  PID: 52688;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 52690;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	5665407	PEPPRO	_RES_

> `Trim_loss_rate`	75.52	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq` (52772)
<pre>
Started analysis of H9_PRO-seq_20_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_20_R1_processed.fastq
Analysis complete for H9_PRO-seq_20_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 0.923GB.  
  PID: 52772;	Command: fastqc;	Return code: 0;	Memory used: 0.179GB

> `FastQC report r1`	fastqc/H9_PRO-seq_20_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_trimmed.fastq`  

> `seqkit rmdup --threads 4 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_noadap.fastq` (52955)
<pre>
[INFO][0m 233822 duplicated records removed
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 0.923GB.  
  PID: 52955;	Command: seqkit;	Return code: 0;	Memory used: 0.499GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_trimmed.fastq` (53108,53109)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.923GB.  
  PID: 53108;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 53109;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	9167377.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	5720857.0	PEPPRO	_RES_

> `Duplicate_reads`	233822.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	49.447	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/processed_R1.flag` (53193)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.923GB.  
  PID: 53193;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed.fastq`  

> `(cutadapt -j 4 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20_R2_cutadapt.txt` (53198)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 0.923GB.  
  PID: 53198;	Command: cutadapt;	Return code: 0;	Memory used: 0.889GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed.fastq` (53413,53414)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 0.923GB.  
  PID: 53413;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 53414;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	11330814	PEPPRO	_RES_

> `Trim_loss_rate`	51.03	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq` (53502)
<pre>
Started analysis of H9_PRO-seq_20_R1_processed.fastq
Approx 5% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 10% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 15% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 20% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 25% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 30% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 35% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 40% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 45% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 50% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 55% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 60% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 65% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 70% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 75% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 80% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 85% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 90% complete for H9_PRO-seq_20_R1_processed.fastq
Approx 95% complete for H9_PRO-seq_20_R1_processed.fastq
Analysis complete for H9_PRO-seq_20_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 0.923GB.  
  PID: 53502;	Command: fastqc;	Return code: 0;	Memory used: 0.175GB

> `FastQC report r1`	fastqc/H9_PRO-seq_20_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed.fastq` (53620)
<pre>
Started analysis of H9_PRO-seq_20_R2_trimmed.fastq
Approx 5% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 10% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 15% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 20% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 25% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 30% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 35% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 40% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 45% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 50% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 55% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 60% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 65% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 70% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 75% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 80% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 85% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 90% complete for H9_PRO-seq_20_R2_trimmed.fastq
Approx 95% complete for H9_PRO-seq_20_R2_trimmed.fastq
Analysis complete for H9_PRO-seq_20_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 0.923GB.  
  PID: 53620;	Command: fastqc;	Return code: 0;	Memory used: 0.182GB

> `FastQC report r2`	fastqc/H9_PRO-seq_20_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20.histogram`  

> `fastq_pair -t 20825402 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_noadap.fastq` (53727)
<pre>
Left paired: 5808026		Right paired: 5808026
Left single: 40785		Right single: 535780
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:45. Running peak memory: 1.008GB.  
  PID: 53727;	Command: fastq_pair;	Return code: 0;	Memory used: 1.008GB


> `flash -q -t 4 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_noadap.fastq.paired.fq -o H9_PRO-seq_20 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt` (53994)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 1.008GB.  
  PID: 53994;	Command: flash;	Return code: 0;	Memory used: 0.058GB


### Plot adapter insertion distribution (06-11 17:59:20) elapsed: 245.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt -u 8` (54314)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 1.008GB.  
  PID: 54314;	Command: Rscript;	Return code: 0;	Memory used: 0.202GB

> `Adapter insertion distribution`	cutadapt/H9_PRO-seq_20_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_PRO-seq_20_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 17:59:27) elapsed: 7.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/cutadapt/H9_PRO-seq_20.hist`

> `Degradation_ratio`	1.0309	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed_dups.fastq` (54377)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 1.008GB.  
  PID: 54377;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/processed_R2.flag` (54414)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.008GB.  
  PID: 54414;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/repaired.flag`  

> `fastq_pair -t 20825402 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed.fastq` (54416)
<pre>
Left paired: 5620420		Right paired: 5620420
Left single: 44987		Right single: 550060
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:42. Running peak memory: 1.056GB.  
  PID: 54416;	Command: fastq_pair;	Return code: 0;	Memory used: 1.056GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq` (55018)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.056GB.  
  PID: 55018;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed.fastq` (55020)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.056GB.  
  PID: 55020;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/repaired.flag` (55024)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.056GB.  
  PID: 55024;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/dups_repaired.flag`  

> `fastq_pair -t 20825402 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed_dups.fastq` (55027)
<pre>
Left paired: 5407657		Right paired: 5407657
Left single: 35362		Right single: 762823
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 1.056GB.  
  PID: 55027;	Command: fastq_pair;	Return code: 0;	Memory used: 1.046GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_trimmed.fastq` (55246)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.056GB.  
  PID: 55246;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed_dups.fastq` (55248)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.056GB.  
  PID: 55248;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/dups_repaired.flag` (55250)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.056GB.  
  PID: 55250;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-11 18:00:54) elapsed: 86.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 18:00:54) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_bt2` (55251)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.056GB.  
  PID: 55251;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq` (55252)
<pre>
</pre>

> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_20 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
5620420 reads; of these:
  5620420 (100.00%) were unpaired; of these:
    5037944 (89.64%) aligned 0 times
    582476 (10.36%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.36% overall alignment rate

> `Aligned_reads_human_rDNA`	1164952.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.28	PEPPRO	_RES_

### Map to human_rDNA (06-11 18:01:39) elapsed: 45.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_dups_bt2` (55654)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.056GB.  
  PID: 55654;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_dups_R2.fq` (55657)
<pre>
</pre>

> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_PRO-seq_20 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
582476 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-11 18:02:35) elapsed: 56.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id H9_PRO-seq_20 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/tmpz5b6rmel -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp.bam` (56007,56009,56010)
<pre>
515625 reads skipped
0 reads lost
5037944 reads; of these:
  5037944 (100.00%) were paired; of these:
    2087736 (41.44%) aligned concordantly 0 times
    2420178 (48.04%) aligned concordantly exactly 1 time
    530030 (10.52%) aligned concordantly >1 times
    ----
    2087736 pairs aligned concordantly 0 times; of these:
      464022 (22.23%) aligned discordantly 1 time
    ----
    1623714 pairs aligned 0 times concordantly or discordantly; of these:
      3247428 mates make up the pairs; of these:
        1333631 (41.07%) aligned 0 times
        699598 (21.54%) aligned exactly 1 time
        1214199 (37.39%) aligned >1 times
86.76% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:22:54. Running peak memory: 3.556GB.  
  PID: 56007;	Command: bowtie2;	Return code: 0;	Memory used: 3.556GB  
  PID: 56009;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 56010;	Command: samtools;	Return code: 0;	Memory used: 0.896GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam` (63957)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 3.556GB.  
  PID: 63957;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	8742257	PEPPRO	_RES_

> `QC_filtered_reads`	5108893	PEPPRO	_RES_

> `Aligned_reads`	3633363.5	PEPPRO	_RES_

> `Alignment_rate`	32.07	PEPPRO	_RES_

> `Total_efficiency`	15.7	PEPPRO	_RES_

> `Read_depth`	1.77	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort_dups.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id H9_PRO-seq_20 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/tmpz5b6rmel -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp_dups.bam` (65706,65707,65708)
<pre>
4892032 reads; of these:
  4892032 (100.00%) were paired; of these:
    1990905 (40.70%) aligned concordantly 0 times
    2383527 (48.72%) aligned concordantly exactly 1 time
    517600 (10.58%) aligned concordantly >1 times
    ----
    1990905 pairs aligned concordantly 0 times; of these:
      457156 (22.96%) aligned discordantly 1 time
    ----
    1533749 pairs aligned 0 times concordantly or discordantly; of these:
      3067498 mates make up the pairs; of these:
        1262056 (41.14%) aligned 0 times
        688075 (22.43%) aligned exactly 1 time
        1117367 (36.43%) aligned >1 times
87.10% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:21:26. Running peak memory: 3.556GB.  
  PID: 65706;	Command: bowtie2;	Return code: 0;	Memory used: 3.544GB  
  PID: 65707;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 65708;	Command: samtools;	Return code: 0;	Memory used: 0.896GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort_dups.bam` (70827)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 3.556GB.  
  PID: 70827;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


### Compress all unmapped read files (06-11 18:52:19) elapsed: 2984.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq` (70864)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.556GB.  
  PID: 70864;	Command: pigz;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R1.fq` (70880)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 3.556GB.  
  PID: 70880;	Command: pigz;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp.bam` (70946)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.556GB.  
  PID: 70946;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	168231	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam` (70961)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.556GB.  
  PID: 70961;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/chr_sizes.bed` (70975,70976,70977,70978)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.556GB.  
  PID: 70976;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 70978;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 70975;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 70977;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/chr_sizes.bed -b -@ 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_noMT.bam` (70981)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.556GB.  
  PID: 70981;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam` (71004)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.556GB.  
  PID: 71004;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam` (71005)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.556GB.  
  PID: 71005;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


### Split BAM file (06-11 18:53:21) elapsed: 62.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam | samtools sort - -@ 4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam` (71011,71012)
<pre>
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 3.556GB.  
  PID: 71011;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 71012;	Command: samtools;	Return code: 0;	Memory used: 1.043GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam | samtools sort - -@ 4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE2.bam` (71157,71158)
<pre>
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 3.556GB.  
  PID: 71157;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 71158;	Command: samtools;	Return code: 0;	Memory used: 0.838GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp_dups.bam` (71219)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.556GB.  
  PID: 71219;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort_dups.bam | samtools sort - -@ 4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_dups_PE1.bam` (71276,71277)
<pre>
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:15) elapsed: 114.0 _TIME_
Child process 55252 (perl) was already terminated.
Child process 55657 (perl) was already terminated.
[W::bgzf_read_block] EOF marker is absent. The input is probably truncated
Child process 71276 (samtools) terminated after 0 sec.
Child process 71277 (samtools) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/recover.lock.aligned_hg38__H9_PRO-seq_20_dups_PE1.bam
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/recover.lock.aligned_hg38__H9_PRO-seq_20_dups_PE2.bam

### Pipeline failed at:  (06-11 18:55:15) elapsed: 0.0 _TIME_

Total time: 1:00:58
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_20 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 4 -M 8000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba26-35c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/
*  Pipeline started at:   (06-11 19:13:13) elapsed: 6.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.8
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16373 insertions(+), 3522 deletions(-)

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `cutadapt`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `4`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE2.fastq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `8000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `H9_PRO-seq_20`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE2.fastq.gz

> `File_mb`	492.8	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:13:14) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1.fastq`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz']

### FASTQ processing:  (06-11 19:13:14) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/processed_R1.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/processed_R2.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/repaired.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/dups_repaired.flag`  

### Prealignments (06-11 19:13:14) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:13:14) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_bt2` (21214)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
  PID: 21214;	Command: mkfifo;	Return code: 0;	Memory used: 0.001GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_bt2

### Map to human_rDNA (06-11 19:13:14) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_dups_bt2` (21215)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.001GB.  
  PID: 21215;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-11 19:13:14) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort_dups.bam`  

### Compress all unmapped read files (06-11 19:13:14) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R1.fq.gz`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq.gz`  

### Split BAM file (06-11 19:13:14) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE2.bam`  

> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort_dups.bam.bai
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/lock.aligned_hg38__H9_PRO-seq_20_dups_PE1.bam
Overwriting target...
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/lock.aligned_hg38__H9_PRO-seq_20_dups_PE2.bam
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort_dups.bam | samtools sort - -@ 4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_dups_PE1.bam` (21221,21222)
<pre>
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 1.05GB.  
  PID: 21221;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 21222;	Command: samtools;	Return code: 0;	Memory used: 1.05GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort_dups.bam | samtools sort - -@ 4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_dups_PE2.bam` (21271,21272)
<pre>
[bam_sort_core] merging from 0 files and 4 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 1.05GB.  
  PID: 21271;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 21272;	Command: samtools;	Return code: 0;	Memory used: 0.846GB


### Calculate library complexity (06-11 19:14:25) elapsed: 71.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_dups_PE1.bam` (21318)
<pre>
BAM_INPUT
TOTAL READS     = 3908887
COUNTS_SUM      = 3908887
DISTINCT READS  = 3.55352e+06
DISTINCT COUNTS = 127
MAX COUNT       = 7844
COUNTS OF 1     = 3.36495e+06
OBSERVED COUNTS (7845)
1	3364952
2	139399
3	26560
4	9493
5	4329
6	2470
7	1478
8	1000
9	680
10	535
11	410
12	334
13	255
14	181
15	166
16	150
17	108
18	101
19	88
20	79
21	68
22	55
23	48
24	47
25	41
26	44
27	36
28	15
29	30
30	19
31	19
32	24
33	17
34	15
35	16
36	11
37	10
38	11
39	13
40	10
41	9
42	7
43	11
44	7
45	10
46	9
47	12
48	4
49	10
50	6
51	2
52	3
53	2
54	2
55	5
56	4
57	2
58	4
59	4
60	5
61	4
62	4
63	2
64	1
65	1
67	2
68	2
70	3
71	2
72	3
73	3
74	3
78	1
79	2
80	1
81	1
82	1
83	1
86	2
89	1
91	1
93	1
94	1
95	1
96	1
98	1
100	1
101	1
104	1
107	1
108	1
109	1
112	1
114	1
115	2
119	1
120	1
125	1
135	1
141	1
143	1
148	1
154	1
155	1
156	1
173	1
192	1
195	1
241	1
251	1
296	1
298	1
303	1
351	1
406	1
414	1
445	1
521	1
529	1
578	1
773	1
1010	1
1984	1
2411	1
2958	1
6314	1
7844	1

sample size: 1000000
sample size: 2000000
sample size: 3000000
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 1.05GB.  
  PID: 21318;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_dups_PE1.bam` (21337)
<pre>
BAM_INPUT
TOTAL READS     = 3908887
DISTINCT READS  = 3.55352e+06
DISTINCT COUNTS = 127
MAX COUNT       = 7844
COUNTS OF 1     = 3.36495e+06
MAX TERMS       = 64
OBSERVED COUNTS (7845)
1	3364952
2	139399
3	26560
4	9493
5	4329
6	2470
7	1478
8	1000
9	680
10	535
11	410
12	334
13	255
14	181
15	166
16	150
17	108
18	101
19	88
20	79
21	68
22	55
23	48
24	47
25	41
26	44
27	36
28	15
29	30
30	19
31	19
32	24
33	17
34	15
35	16
36	11
37	10
38	11
39	13
40	10
41	9
42	7
43	11
44	7
45	10
46	9
47	12
48	4
49	10
50	6
51	2
52	3
53	2
54	2
55	5
56	4
57	2
58	4
59	4
60	5
61	4
62	4
63	2
64	1
65	1
67	2
68	2
70	3
71	2
72	3
73	3
74	3
78	1
79	2
80	1
81	1
82	1
83	1
86	2
89	1
91	1
93	1
94	1
95	1
96	1
98	1
100	1
101	1
104	1
107	1
108	1
109	1
112	1
114	1
115	2
119	1
120	1
125	1
135	1
141	1
143	1
148	1
154	1
155	1
156	1
173	1
192	1
195	1
241	1
251	1
296	1
298	1
303	1
351	1
406	1
414	1
445	1
521	1
529	1
578	1
773	1
1010	1
1984	1
2411	1
2958	1
6314	1
7844	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
_................................_....................................................................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 1.05GB.  
  PID: 21337;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_counts.txt` (21564)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 1.05GB.  
  PID: 21564;	Command: echo;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_plot` (21575)
<pre>
Processing H9_PRO-seq_20
INFO: Found real counts for H9_PRO-seq_20 - Total (M): 3.896875 Unique (M): 3.908887

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 1.05GB.  
  PID: 21575;	Command: Rscript;	Return code: 0;	Memory used: 0.129GB

> `Library complexity`	QC_hg38/H9_PRO-seq_20_preseq_plot.pdf	Library complexity	QC_hg38/H9_PRO-seq_20_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8521	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-11 19:15:30) elapsed: 64.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam` (21598)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 1.05GB.  
  PID: 21598;	Command: samtools;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam -c 4 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_bamQC.tsv` (21602)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/tmp_H9_PRO-seq_20_PE1_qmolbs_6'
Processing with 4 cores...
Discarding 114 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270725v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 81 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 1.05GB.  
  PID: 21602;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.285GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	1948437.5	PEPPRO	_RES_

> `PBC2`	1948437.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_unmap.bam`  

> `samtools view -b -@ 4 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_unmap.bam` (21625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 1.05GB.  
  PID: 21625;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_temp.bam`

> `Unmapped_reads`	1333631	PEPPRO	_RES_

### Split BAM by strand (06-11 19:15:51) elapsed: 21.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam` (21644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 1.05GB.  
  PID: 21644;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam` (21662)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 1.05GB.  
  PID: 21662;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-11 19:16:24) elapsed: 32.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (21677)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.05GB.  
  PID: 21677;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_plus_TssEnrichment.txt` (21678)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 1.05GB.  
  PID: 21678;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.321GB


> `TSS_coding_score`	33.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_minus_TssEnrichment.txt` (21694)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 1.05GB.  
  PID: 21694;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.354GB


> `TSS_non-coding_score`	12.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_minus_TssEnrichment.txt` (21710)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 1.05GB.  
  PID: 21710;	Command: Rscript;	Return code: 0;	Memory used: 0.285GB

> `TSS enrichment`	QC_hg38/H9_PRO-seq_20_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_20_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt` (21728,21729,21730,21731)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.05GB.  
  PID: 21728;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 21730;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 21729;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 21731;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_keep.txt` (21733)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.05GB.  
  PID: 21733;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-11 19:16:42) elapsed: 18.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_ensembl_tss.bed` (21735,21736)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 1.05GB.  
  PID: 21735;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 21736;	Command: bedtools;	Return code: 0;	Memory used: 0.099GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_ensembl_gene_body.bed` (21740,21741)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.05GB.  
  PID: 21740;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 21741;	Command: bedtools;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_TSS_density.bed` (21743,21744,21745,21746)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 1.05GB.  
  PID: 21744;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 21746;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 21743;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 21745;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_gene_body_density.bed` (21754,21755,21756)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 1.05GB.  
  PID: 21755;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 21754;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB  
  PID: 21756;	Command: sort;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/tmpnoh0qgsn` (21768,21769,21770)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 1.05GB.  
  PID: 21768;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 21770;	Command: env;	Return code: 0;	Memory used: 0.005GB  
  PID: 21769;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/tmpnoh0qgsn | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-11 19:17:00) elapsed: 18.0 _TIME_

Total time: 0:03:53
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/tmpnoh0qgsn | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_PRO-seq_20 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 4 -M 8000 --input2 /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-aj37-17c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/
*  Pipeline started at:   (06-14 21:18:39) elapsed: 5.0 _TIME_

### Version log:

*       Python version:  3.6.6
*          Pypiper dir:  `/sfs/qumulo/qhome/jps3dp/.local/lib/python3.6/site-packages/pypiper`
*      Pypiper version:  0.12.1
*         Pipeline dir:  `/sfs/lustre/bahamut/scratch/jps3dp/tools/databio/peppro/pipelines`
*     Pipeline version:  0.9.8
*        Pipeline hash:  887a9a85408962dc3ca070dc954e59a5d4d73a10
*      Pipeline branch:  * master
*        Pipeline date:  2020-06-09 11:36:49 -0400
*        Pipeline diff:  40 files changed, 16373 insertions(+), 3522 deletions(-)

### Arguments passed to pipeline:

*           `TSS_name`:  `None`
*            `adapter`:  `cutadapt`
*          `anno_name`:  `None`
*         `complexity`:  `False`
*        `config_file`:  `peppro.yaml`
*              `cores`:  `4`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE2.fastq.gz']`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `8000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `True`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `H9_PRO-seq_20`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_PRO-seq_20pct_PE2.fastq.gz

> `File_mb`	492.8	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:18:40) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R1.fastq`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/H9_PRO-seq_20_R2.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/H9_PRO-seq_20_R2.fastq.gz']

### FASTQ processing:  (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/processed_R1.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/processed_R2.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/repaired.flag`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/fastq/dups_repaired.flag`  

### Prealignments (06-14 21:18:40) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:18:40) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_bt2` (315045)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 315045;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_bt2

### Map to human_rDNA (06-14 21:18:40) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/human_rDNA_dups_bt2` (315138)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 315138;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq.gz`  
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_sort.bam`  

### Compress all unmapped read files (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R2.fq.gz`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/prealignments/H9_PRO-seq_20_human_rDNA_unmap_R1.fq.gz`  

### Split BAM file (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE2.bam`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_unmap.bam`  

### Split BAM by strand (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam`  

### Calculate TSS enrichment (06-14 21:18:40) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/H9_PRO-seq_20_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_PRO-seq_20_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:18:40) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/tmp776j47gc` (315371,315391,315420)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 315371;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 315420;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 315391;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/tmp776j47gc | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0047564) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/tmp776j47gc > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_pause_index.bed` (315958)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 315958;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	33.2	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_pause_index.bed` (316144)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.204GB.  
  PID: 316144;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Pause index`	QC_hg38/H9_PRO-seq_20_pause_index.pdf	Pause index	QC_hg38/H9_PRO-seq_20_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_pause_index.bed` (325365)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 325365;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:18:47) elapsed: 7.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam`
3633363.5 1411171

> `Plus_FRiP`	0.39	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam`
3633363.5 1313835

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_gene_sort.bed` (331058,331060)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331058;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 331060;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_gene_coverage.bed` (331425)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 0.204GB.  
  PID: 331425;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/hg38_annotations.bed.gz` (331791)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331791;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/hg38_annotations.bed` (331792)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331792;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:19:04) elapsed: 16.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/raw/hg38_annotations.bed` (331800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331800;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Enhancer_sort.bed` (331810,331811,331812,331813)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331810;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331811;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 331813;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 331812;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Enhancer_plus_coverage.bed` (331818)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331818;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Enhancer_minus_coverage.bed` (331830)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331830;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_sort.bed` (331846,331847,331848,331849)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331846;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331847;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 331849;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 331848;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_plus_coverage.bed` (331851)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331851;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_minus_coverage.bed` (331865)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331865;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_Flanking_Region"` (331878)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331878;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed` (331879,331880,331881,331882)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331879;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331881;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 331880;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 331882;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_Flanking_Region_plus_coverage.bed` (331886)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331886;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_Flanking_Region_minus_coverage.bed` (331896)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331896;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/5_UTR"` (331902)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331902;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/5_UTR_sort.bed` (331903,331904,331905,331906)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331903;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331904;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 331906;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 331905;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_5_UTR_plus_coverage.bed` (331912)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331912;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_5_UTR_minus_coverage.bed` (331923)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331923;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/3_UTR"` (331931)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 331931;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/3_UTR_sort.bed` (331932,331933,331934,331935)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331932;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331933;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 331935;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB  
  PID: 331934;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_3_UTR_plus_coverage.bed` (331939)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331939;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_3_UTR_minus_coverage.bed` (331949)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331949;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Exon_sort.bed` (331971,331972,331973,331974)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.204GB.  
  PID: 331971;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331972;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 331974;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 331973;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Exon_plus_coverage.bed` (331979)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 331979;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Exon_minus_coverage.bed` (331984)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 331984;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Intron_sort.bed` (331989,331990,331991,331992)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 331989;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 331991;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 331990;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 331992;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Intron_plus_coverage.bed` (331995)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.204GB.  
  PID: 331995;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Intron_minus_coverage.bed` (332002)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.204GB.  
  PID: 332002;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB


### Plot cFRiF/FRiF (06-14 21:19:58) elapsed: 55.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_20 -z 3099922541 -n 1990532 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Intron_plus_coverage.bed` (332013)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 0.443GB.  
  PID: 332013;	Command: Rscript;	Return code: 0;	Memory used: 0.443GB

> `cFRiF`	QC_hg38/H9_PRO-seq_20_cFRiF.pdf	cFRiF	QC_hg38/H9_PRO-seq_20_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_PRO-seq_20 -z 3099922541 -n 1990532 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_Intron_plus_coverage.bed` (332332)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 0.443GB.  
  PID: 332332;	Command: Rscript;	Return code: 0;	Memory used: 0.443GB

> `FRiF`	QC_hg38/H9_PRO-seq_20_FRiF.pdf	FRiF	QC_hg38/H9_PRO-seq_20_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 21:21:02) elapsed: 64.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_exons_sort.bed` (332406,332407)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.443GB.  
  PID: 332407;	Command: bedtools;	Return code: 0;	Memory used: 0.087GB  
  PID: 332406;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_introns_sort.bed` (332428,332429,332430)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.443GB.  
  PID: 332428;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 332430;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 332429;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exons_coverage.bed` (332448)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.443GB.  
  PID: 332448;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_introns_coverage.bed` (332455)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.443GB.  
  PID: 332455;	Command: bedtools;	Return code: 0;	Memory used: 0.029GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/3.6333635)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exons_rpkm.bed` (332462,332463,332464)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.443GB.  
  PID: 332462;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 332464;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 332463;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/3.6333635)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_introns_rpkm.bed` (332467,332468,332469)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.443GB.  
  PID: 332467;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 332469;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 332468;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exon_intron_ratios.bed` (332472,332473,332474)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.443GB.  
  PID: 332472;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 332474;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 332473;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.31	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exon_intron_ratios.bed --annotate` (332480)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.443GB.  
  PID: 332480;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/H9_PRO-seq_20_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_PRO-seq_20_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/QC_hg38/H9_PRO-seq_20_exon_intron_ratios.bed` (332501)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.443GB.  
  PID: 332501;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (06-14 21:21:31) elapsed: 29.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam` (332509)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.443GB.  
  PID: 332509;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge --scale 3633363.5` (332515)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_plus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_20_plus_cuttrace_el5rfk2l'
Processing with 1 cores...
Discarding 132 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270714v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr14_GL000009v2_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr15_KI270727v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000218v1']
Keeping 63 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270726v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270539v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrEBV']
Reduce step (merge files)...
Merging 63 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_plus_exact_body_0-mer.bw'
Merging 63 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:16:49. Running peak memory: 1.33GB.  
  PID: 332515;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.33GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam` (135538)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 1.33GB.  
  PID: 135538;	Command: samtools;	Return code: 0;	Memory used: 0.007GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_minus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge --scale 3633363.5` (135541)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/aligned_hg38/H9_PRO-seq_20_minus.bam'
Temporary files will be stored in: 'tmp_H9_PRO-seq_20_minus_cuttrace_g_w9nf7w'
Processing with 1 cores...
stdin is empty of data
Discarding 131 chunk(s) of reads: ['chrM', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_KI270722v1_random', 'chr14_KI270723v1_random', 'chr14_KI270725v1_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1']
Keeping 64 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr14_GL000225v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270737v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 64 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_minus_exact_body_0-mer.bw'
Merging 64 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_PRO-seq_20/signal_hg38/H9_PRO-seq_20_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:16:35. Running peak memory: 1.33GB.  
  PID: 135541;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.281GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  0:36:25
*  Total elapsed time (all runs):  4:44:32
*         Peak memory (this run):  1.3301 GB
*        Pipeline completed time: 2020-06-14 21:54:59
