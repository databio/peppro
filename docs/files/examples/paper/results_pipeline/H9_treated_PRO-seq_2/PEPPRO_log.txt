### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_treated_PRO-seq_2 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep2_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep2_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-28c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/
*  Pipeline started at:   (06-15 07:17:18) elapsed: 1.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep2_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep2_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_treated_PRO-seq_2`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep2_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep2_PE2.fastq.gz

> `File_mb`	2737.32	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:19) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep2_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R1.fastq.gz` (138709)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 138709;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep2_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R2.fastq.gz` (138712)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 138712;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1.fastq` (138715)
<pre>
</pre>
Command completed. Elapsed time: 0:02:49. Running peak memory: 0.002GB.  
  PID: 138715;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2.fastq` (45196)
<pre>
</pre>
Command completed. Elapsed time: 0:01:47. Running peak memory: 0.002GB.  
  PID: 45196;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	114083518	PEPPRO	_RES_

> `Fastq_reads`	114083518	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/H9_treated_PRO-seq_2_R2.fastq.gz']

### FASTQ processing:  (06-15 07:25:43) elapsed: 504.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R1_cutadapt.txt` (3602)
<pre>
</pre>
Command completed. Elapsed time: 0:02:40. Running peak memory: 3.232GB.  
  PID: 3602;	Command: cutadapt;	Return code: 0;	Memory used: 3.232GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq` (3786,3787)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 3.232GB.  
  PID: 3786;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 3787;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	28674037	PEPPRO	_RES_

> `Trim_loss_rate`	74.87	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq` (3848)
<pre>
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
Command completed. Elapsed time: 0:00:58. Running peak memory: 3.232GB.  
  PID: 3848;	Command: fastqc;	Return code: 0;	Memory used: 0.182GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_2_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq` (3943)
<pre>
[INFO][0m 2388789 duplicated records removed
</pre>
Command completed. Elapsed time: 0:01:04. Running peak memory: 3.232GB.  
  PID: 3943;	Command: seqkit;	Return code: 0;	Memory used: 2.036GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq` (4301,4302)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 3.232GB.  
  PID: 4301;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 4302;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	43984199.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	27614116.0	PEPPRO	_RES_

> `Duplicate_reads`	2388789.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	48.4104	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/processed_R1.flag` (4515)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.232GB.  
  PID: 4515;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_R2_cutadapt.txt` (4517)
<pre>
</pre>
Command completed. Elapsed time: 0:02:10. Running peak memory: 3.35GB.  
  PID: 4517;	Command: cutadapt;	Return code: 0;	Memory used: 3.35GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq` (4928,4929)
<pre>
</pre>
Command completed. Elapsed time: 0:00:47. Running peak memory: 3.35GB.  
  PID: 4928;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 4929;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	57348074	PEPPRO	_RES_

> `Trim_loss_rate`	49.73	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq` (5021)
<pre>
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
Command completed. Elapsed time: 0:00:59. Running peak memory: 3.35GB.  
  PID: 5021;	Command: fastqc;	Return code: 0;	Memory used: 0.179GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_2_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq` (5107)
<pre>
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
Command completed. Elapsed time: 0:01:03. Running peak memory: 3.35GB.  
  PID: 5107;	Command: fastqc;	Return code: 0;	Memory used: 0.17GB

> `FastQC report r2`	fastqc/H9_treated_PRO-seq_2_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.histogram`  

> `fastq_pair -t 102675166 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq` (5194)
<pre>
Left paired: 29252476		Right paired: 29252476
Left single: 175167		Right single: 2238349
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:04:07. Running peak memory: 6.557GB.  
  PID: 5194;	Command: fastq_pair;	Return code: 0;	Memory used: 6.557GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_noadap.fastq.paired.fq -o H9_treated_PRO-seq_2 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt` (5825)
<pre>
</pre>
Command completed. Elapsed time: 0:01:04. Running peak memory: 6.557GB.  
  PID: 5825;	Command: flash;	Return code: 0;	Memory used: 0.117GB


### Plot adapter insertion distribution (06-15 07:44:09) elapsed: 1106.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt -u 8` (6026)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.557GB.  
  PID: 6026;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/H9_treated_PRO-seq_2_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_treated_PRO-seq_2_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:44:15) elapsed: 6.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/cutadapt/H9_treated_PRO-seq_2.hist`

> `Degradation_ratio`	1.0037	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq` (6203)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 6.557GB.  
  PID: 6203;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/processed_R2.flag` (6234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 6234;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/repaired.flag`  

> `fastq_pair -t 102675166 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq` (6235)
<pre>
Left paired: 28465823		Right paired: 28465823
Left single: 208214		Right single: 2295865
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:40. Running peak memory: 6.557GB.  
  PID: 6235;	Command: fastq_pair;	Return code: 0;	Memory used: 5.268GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq` (6592)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 6592;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq` (6595)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 6595;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/repaired.flag` (6596)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 6596;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/dups_repaired.flag`  

> `fastq_pair -t 102675166 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq` (6597)
<pre>
Left paired: 26273211		Right paired: 26273211
Left single: 141535		Right single: 4488477
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:34. Running peak memory: 6.557GB.  
  PID: 6597;	Command: fastq_pair;	Return code: 0;	Memory used: 5.269GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq` (7001)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 7001;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq` (7003)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 7003;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/dups_repaired.flag` (7004)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 7004;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:50:05) elapsed: 350.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:50:05) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_bt2` (7005)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 7005;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq` (7006)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_2 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
28465823 reads; of these:
  28465823 (100.00%) were unpaired; of these:
    25736513 (90.41%) aligned 0 times
    2729310 (9.59%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.59% overall alignment rate

> `Aligned_reads_human_rDNA`	5458620.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.52	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:54:07) elapsed: 242.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_dups_bt2` (7278)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 7278;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_dups_R2.fq` (7279)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_2 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/fastq/H9_treated_PRO-seq_2_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
2729310 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:57:59) elapsed: 231.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_2 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/tmpj5bu8aqe -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam` (7724,7725,7726)
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
[bam_sort_core] merging from 15 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:49:01. Running peak memory: 6.557GB.  
  PID: 7724;	Command: bowtie2;	Return code: 0;	Memory used: 3.778GB  
  PID: 7725;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 7726;	Command: samtools;	Return code: 0;	Memory used: 0.906GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam` (13046)
<pre>
</pre>
Command completed. Elapsed time: 0:01:33. Running peak memory: 6.557GB.  
  PID: 13046;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	44741250	PEPPRO	_RES_

> `QC_filtered_reads`	25860518	PEPPRO	_RES_

> `Aligned_reads`	18880731.5	PEPPRO	_RES_

> `Alignment_rate`	32.92	PEPPRO	_RES_

> `Total_efficiency`	16.55	PEPPRO	_RES_

> `Read_depth`	3.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_2 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/tmpj5bu8aqe -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp_dups.bam` (14640,14642,14643)
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
[bam_sort_core] merging from 14 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:44:45. Running peak memory: 6.557GB.  
  PID: 14640;	Command: bowtie2;	Return code: 0;	Memory used: 3.765GB  
  PID: 14642;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 14643;	Command: samtools;	Return code: 0;	Memory used: 0.905GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam` (19264)
<pre>
</pre>
Command completed. Elapsed time: 0:01:23. Running peak memory: 6.557GB.  
  PID: 19264;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 09:47:18) elapsed: 6560.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R1.fq` (19418)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 6.557GB.  
  PID: 19418;	Command: pigz;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/prealignments/H9_treated_PRO-seq_2_human_rDNA_unmap_R2.fq` (19451)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 6.557GB.  
  PID: 19451;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam` (19488)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 6.557GB.  
  PID: 19488;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	670904	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam` (19546)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.557GB.  
  PID: 19546;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/chr_sizes.bed` (19574,19575,19576,19577)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 19575;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 19577;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 19574;	Command: samtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 19576;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_noMT.bam` (19579)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 6.557GB.  
  PID: 19579;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam` (19690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 19690;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam` (19691)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.557GB.  
  PID: 19691;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


### Split BAM file (06-15 09:50:28) elapsed: 190.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam` (19978,19979)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:51. Running peak memory: 6.557GB.  
  PID: 19978;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 19979;	Command: samtools;	Return code: 0;	Memory used: 5.526GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE2.bam` (20253,20254)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:21. Running peak memory: 6.557GB.  
  PID: 20253;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20254;	Command: samtools;	Return code: 0;	Memory used: 4.484GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp_dups.bam` (20756)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 6.557GB.  
  PID: 20756;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE1.bam` (20791,20792)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:46. Running peak memory: 6.557GB.  
  PID: 20791;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20792;	Command: samtools;	Return code: 0;	Memory used: 5.349GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE2.bam` (21036,21037)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:11. Running peak memory: 6.557GB.  
  PID: 21036;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 21037;	Command: samtools;	Return code: 0;	Memory used: 4.334GB


### Calculate library complexity (06-15 10:02:11) elapsed: 703.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE1.bam` (21548)
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
  PID: 21548;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE1.bam` (21677)
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
..._......................................................_.................._.........................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:53. Running peak memory: 6.557GB.  
  PID: 21677;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_counts.txt` (22125)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.557GB.  
  PID: 22125;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_plot` (22230)
<pre>
Processing H9_treated_PRO-seq_2
INFO: Found real counts for H9_treated_PRO-seq_2 - Total (M): 20.247016 Unique (M): 19.595024

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.557GB.  
  PID: 22230;	Command: Rscript;	Return code: 0;	Memory used: 0.285GB

> `Library complexity`	QC_hg38/H9_treated_PRO-seq_2_preseq_plot.pdf	Library complexity	QC_hg38/H9_treated_PRO-seq_2_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8189	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 10:06:27) elapsed: 256.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam` (22265)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.557GB.  
  PID: 22265;	Command: samtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_bamQC.tsv` (22287)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/tmp_H9_treated_PRO-seq_2_PE1_0_3why97'
Processing with 12 cores...
Discarding 98 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270757v1']
Keeping 97 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270372v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 6.557GB.  
  PID: 22287;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.242GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	10123508.0	PEPPRO	_RES_

> `PBC2`	10123508.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_unmap.bam` (22340)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 6.557GB.  
  PID: 22340;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_temp.bam`

> `Unmapped_reads`	6731776	PEPPRO	_RES_

### Split BAM by strand (06-15 10:07:22) elapsed: 55.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam` (22402)
<pre>
</pre>
Command completed. Elapsed time: 0:01:08. Running peak memory: 6.557GB.  
  PID: 22402;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam` (22482)
<pre>
</pre>
Command completed. Elapsed time: 0:01:07. Running peak memory: 6.557GB.  
  PID: 22482;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 10:09:36) elapsed: 134.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (22655)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 22655;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_plus_TssEnrichment.txt` (22656)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.557GB.  
  PID: 22656;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.884GB


> `TSS_coding_score`	57.4	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_minus_TssEnrichment.txt` (22689)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.557GB.  
  PID: 22689;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.887GB


> `TSS_non-coding_score`	16.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_minus_TssEnrichment.txt` (22724)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.557GB.  
  PID: 22724;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/H9_treated_PRO-seq_2_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_treated_PRO-seq_2_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt` (22746,22747,22748,22749)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 22746;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 22748;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 22747;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 22749;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt` (22751)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 22751;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 10:09:57) elapsed: 21.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed` (22753,22754)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.557GB.  
  PID: 22753;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 22754;	Command: bedtools;	Return code: 0;	Memory used: 0.097GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed` (22757,22758)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 22757;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 22758;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_TSS_density.bed` (22761,22762,22763,22764)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.557GB.  
  PID: 22762;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 22764;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 22761;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 22763;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_gene_body_density.bed` (23067,23068,23069)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.557GB.  
  PID: 23068;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 23067;	Command: bedtools;	Return code: 0;	Memory used: 0.058GB  
  PID: 23069;	Command: sort;	Return code: 0;	Memory used: 0.003GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/tmpf0akdhol` (23163,23164,23165)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 23163;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 23165;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 23164;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/tmpf0akdhol | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0199168) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/tmpf0akdhol > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed` (23172)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 23172;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	59.52	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed` (23177)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.557GB.  
  PID: 23177;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/H9_treated_PRO-seq_2_pause_index.pdf	Pause index	QC_hg38/H9_treated_PRO-seq_2_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_pause_index.bed` (23204)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 23204;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 10:11:02) elapsed: 66.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam`
18880731.5 7474169

> `Plus_FRiP`	0.4	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam`
18880731.5 7075799

> `Minus_FRiP`	0.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_gene_sort.bed` (23276,23277)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 23276;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 23277;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_gene_coverage.bed` (23279)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 6.557GB.  
  PID: 23279;	Command: bedtools;	Return code: 0;	Memory used: 0.061GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed.gz` (23328)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 23328;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed` (23329)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 23329;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 10:12:09) elapsed: 67.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/raw/hg38_annotations.bed` (23338)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 23338;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer_sort.bed` (23340,23341,23342,23343)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 23340;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23341;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 23343;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 23342;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_plus_coverage.bed` (23345)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 23345;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_minus_coverage.bed` (23359)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.557GB.  
  PID: 23359;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_sort.bed` (23405,23406,23407,23408)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 23405;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23406;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 23408;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 23407;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_plus_coverage.bed` (23410)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.557GB.  
  PID: 23410;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_minus_coverage.bed` (23498)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.557GB.  
  PID: 23498;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region"` (23519)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 23519;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed` (23520,23521,23522,23523)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 23520;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23522;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 23521;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 23523;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed` (23526)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 23526;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_minus_coverage.bed` (23551)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 23551;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR"` (23590)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 23590;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR_sort.bed` (23591,23592,23593,23594)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 23591;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23592;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 23594;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 23593;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_plus_coverage.bed` (23597)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 23597;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_minus_coverage.bed` (23629)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 23629;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR"` (23650)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 23650;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR_sort.bed` (23652,23653,23654,23655)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 23652;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23653;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 23655;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 23654;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_plus_coverage.bed` (23666)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 23666;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_minus_coverage.bed` (23696)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.557GB.  
  PID: 23696;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon_sort.bed` (23720,23721,23722,23723)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.557GB.  
  PID: 23720;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23721;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 23723;	Command: bedtools;	Return code: 0;	Memory used: 0.169GB  
  PID: 23722;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_plus_coverage.bed` (23736)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.557GB.  
  PID: 23736;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_minus_coverage.bed` (23758)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.557GB.  
  PID: 23758;	Command: bedtools;	Return code: 0;	Memory used: 0.017GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron_sort.bed` (23791,23792,23793,23794)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 23791;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23793;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 23792;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 23794;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_plus_coverage.bed` (23797)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 6.557GB.  
  PID: 23797;	Command: bedtools;	Return code: 0;	Memory used: 0.081GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_minus_coverage.bed` (24580)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.557GB.  
  PID: 24580;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB


### Plot cFRiF/FRiF (06-15 10:15:29) elapsed: 200.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_2 -z 3099922541 -n 10264898 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_plus_coverage.bed` (24620)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 6.557GB.  
  PID: 24620;	Command: Rscript;	Return code: 0;	Memory used: 0.443GB

> `cFRiF`	QC_hg38/H9_treated_PRO-seq_2_cFRiF.pdf	cFRiF	QC_hg38/H9_treated_PRO-seq_2_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_2 -z 3099922541 -n 10264898 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_Intron_plus_coverage.bed` (24706)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.557GB.  
  PID: 24706;	Command: Rscript;	Return code: 0;	Memory used: 0.451GB

> `FRiF`	QC_hg38/H9_treated_PRO-seq_2_FRiF.pdf	FRiF	QC_hg38/H9_treated_PRO-seq_2_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 10:16:34) elapsed: 65.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_exons_sort.bed` (24751,24752)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.557GB.  
  PID: 24752;	Command: bedtools;	Return code: 0;	Memory used: 0.094GB  
  PID: 24751;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_introns_sort.bed` (24776,24777,24779)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.557GB.  
  PID: 24776;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 24779;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 24777;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_coverage.bed` (24789)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.557GB.  
  PID: 24789;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_coverage.bed` (24819)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.557GB.  
  PID: 24819;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.8807315)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_rpkm.bed` (24855,24856,24857)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 24855;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 24857;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 24856;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.8807315)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_rpkm.bed` (24859,24860,24861)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.557GB.  
  PID: 24859;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 24861;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 24860;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed` (24864,24865,24866)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 24864;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 24866;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 24865;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.16	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed --annotate` (24872)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.557GB.  
  PID: 24872;	Command: Rscript;	Return code: 0;	Memory used: 0.302GB

> `mRNA contamination`	QC_hg38/H9_treated_PRO-seq_2_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_treated_PRO-seq_2_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/QC_hg38/H9_treated_PRO-seq_2_exon_intron_ratios.bed` (24893)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.557GB.  
  PID: 24893;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 10:17:45) elapsed: 71.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam` (24901)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.557GB.  
  PID: 24901;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 18880731.5` (24908)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_plus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_2_plus_cuttrace_yuuf674x'
Processing with 4 cores...
stdin is empty of data
Discarding 111 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr5_GL000208v1_random', 'chr9_KI270720v1_random', 'chr14_KI270723v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 84 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr4_GL000008v2_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270468v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270372v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_exact_body_0-mer.bw'
Merging 84 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:24. Running peak memory: 6.557GB.  
  PID: 24908;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.724GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam` (28151)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.557GB.  
  PID: 28151;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 18880731.5` (28158)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/aligned_hg38/H9_treated_PRO-seq_2_minus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_2_minus_cuttrace_w7kmqi2c'
Processing with 4 cores...
Discarding 103 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr22_KI270736v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270757v1']
Keeping 92 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 92 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_exact_body_0-mer.bw'
Merging 92 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_2/signal_hg38/H9_treated_PRO-seq_2_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:16. Running peak memory: 6.557GB.  
  PID: 28158;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.65GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:15:23
*  Total elapsed time (all runs):  5:47:25
*         Peak memory (this run):  6.5575 GB
*        Pipeline completed time: 2020-06-15 10:32:41
