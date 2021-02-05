### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_treated_PRO-seq_1 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-28c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/
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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_treated_PRO-seq_1`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz

> `File_mb`	2710.39	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:19) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep1_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz` (138708)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 138708;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep1_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz` (138711)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 138711;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1.fastq` (138714)
<pre>
</pre>
Command completed. Elapsed time: 0:02:34. Running peak memory: 0.002GB.  
  PID: 138714;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2.fastq` (31182)
<pre>
</pre>
Command completed. Elapsed time: 0:01:45. Running peak memory: 0.002GB.  
  PID: 31182;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	111300168	PEPPRO	_RES_

> `Fastq_reads`	111300168	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/H9_treated_PRO-seq_1_R2.fastq.gz']

### FASTQ processing:  (06-15 07:24:45) elapsed: 446.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt` (3305)
<pre>
</pre>
Command completed. Elapsed time: 0:01:45. Running peak memory: 3.069GB.  
  PID: 3305;	Command: cutadapt;	Return code: 0;	Memory used: 3.069GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (3668,3669)
<pre>
</pre>
Command completed. Elapsed time: 0:01:16. Running peak memory: 3.069GB.  
  PID: 3669;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB  
  PID: 3668;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB

Evaluating read trimming

> `Trimmed_reads`	29119715	PEPPRO	_RES_

> `Trim_loss_rate`	73.84	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (3741)
<pre>
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
Command completed. Elapsed time: 0:01:02. Running peak memory: 3.069GB.  
  PID: 3741;	Command: fastqc;	Return code: 0;	Memory used: 0.179GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq` (3820)
<pre>
[INFO][0m 2891628 duplicated records removed
</pre>
Command completed. Elapsed time: 0:01:12. Running peak memory: 3.069GB.  
  PID: 3820;	Command: seqkit;	Return code: 0;	Memory used: 1.444GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq` (4205,4206)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 3.069GB.  
  PID: 4205;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 4206;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	43631018.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	25201427.0	PEPPRO	_RES_

> `Duplicate_reads`	2891628.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	45.2855	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/processed_R1.flag` (4363)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.069GB.  
  PID: 4363;	Command: touch;	Return code: 0;	Memory used: 0.0GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_R2_cutadapt.txt` (4366)
<pre>
</pre>
Command completed. Elapsed time: 0:01:41. Running peak memory: 3.14GB.  
  PID: 4366;	Command: cutadapt;	Return code: 0;	Memory used: 3.14GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (4511,4512)
<pre>
</pre>
Command completed. Elapsed time: 0:01:04. Running peak memory: 3.14GB.  
  PID: 4512;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB  
  PID: 4511;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB

Evaluating read trimming

> `Trimmed_reads`	58239430	PEPPRO	_RES_

> `Trim_loss_rate`	47.67	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (4880)
<pre>
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
Command completed. Elapsed time: 0:01:04. Running peak memory: 3.14GB.  
  PID: 4880;	Command: fastqc;	Return code: 0;	Memory used: 0.176GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_1_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (4962)
<pre>
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
Command completed. Elapsed time: 0:01:22. Running peak memory: 3.14GB.  
  PID: 4962;	Command: fastqc;	Return code: 0;	Memory used: 0.162GB

> `FastQC report r2`	fastqc/H9_treated_PRO-seq_1_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.histogram`  

> `fastq_pair -t 100170151 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq` (5090)
<pre>
Left paired: 30307647		Right paired: 30307647
Left single: 141010		Right single: 2158148
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:03:15. Running peak memory: 6.727GB.  
  PID: 5090;	Command: fastq_pair;	Return code: 0;	Memory used: 6.727GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_noadap.fastq.paired.fq -o H9_treated_PRO-seq_1 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt` (5534)
<pre>
</pre>
Command completed. Elapsed time: 0:01:26. Running peak memory: 6.727GB.  
  PID: 5534;	Command: flash;	Return code: 0;	Memory used: 0.115GB


### Plot adapter insertion distribution (06-15 07:42:16) elapsed: 1051.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt -u 8` (5753)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.727GB.  
  PID: 5753;	Command: Rscript;	Return code: 0;	Memory used: 0.112GB

> `Adapter insertion distribution`	cutadapt/H9_treated_PRO-seq_1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_treated_PRO-seq_1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:42:24) elapsed: 8.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/cutadapt/H9_treated_PRO-seq_1.hist`

> `Degradation_ratio`	1.158	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq` (5784)
<pre>
</pre>
Command completed. Elapsed time: 0:00:23. Running peak memory: 6.727GB.  
  PID: 5784;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/processed_R2.flag` (5807)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 5807;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/repaired.flag`  

> `fastq_pair -t 100170151 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (5809)
<pre>
Left paired: 28948010		Right paired: 28948010
Left single: 171705		Right single: 2245320
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:03:04. Running peak memory: 6.727GB.  
  PID: 5809;	Command: fastq_pair;	Return code: 0;	Memory used: 6.48GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq` (6506)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 6506;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq` (6508)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 6508;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/repaired.flag` (6509)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 6509;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/dups_repaired.flag`  

> `fastq_pair -t 100170151 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq` (6511)
<pre>
Left paired: 26293117		Right paired: 26293117
Left single: 111289		Right single: 4900213
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:27. Running peak memory: 6.727GB.  
  PID: 6511;	Command: fastq_pair;	Return code: 0;	Memory used: 5.33GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq` (6640)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 6640;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq` (6642)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 6642;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/dups_repaired.flag` (6643)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 6643;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:48:21) elapsed: 357.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:48:21) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2` (6644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 6644;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq` (6645)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
28948010 reads; of these:
  28948010 (100.00%) were unpaired; of these:
    25898407 (89.47%) aligned 0 times
    3049603 (10.53%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.53% overall alignment rate

> `Aligned_reads_human_rDNA`	6099206.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.47	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:51:54) elapsed: 213.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2` (7114)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 7114;	Command: mkfifo;	Return code: 0;	Memory used: 0.002GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R2.fq` (7115)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_1 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/fastq/H9_treated_PRO-seq_1_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
3049603 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:55:12) elapsed: 198.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_1 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/tmpc6_9h7v5 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam` (7546,7547,7548)
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
[bam_sort_core] merging from 15 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:49:39. Running peak memory: 6.727GB.  
  PID: 7547;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 7546;	Command: bowtie2;	Return code: 0;	Memory used: 3.798GB  
  PID: 7548;	Command: samtools;	Return code: 0;	Memory used: 0.897GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (12725)
<pre>
</pre>
Command completed. Elapsed time: 0:01:39. Running peak memory: 6.727GB.  
  PID: 12725;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	44081295	PEPPRO	_RES_

> `QC_filtered_reads`	26072183	PEPPRO	_RES_

> `Aligned_reads`	18009111.5	PEPPRO	_RES_

> `Alignment_rate`	30.92	PEPPRO	_RES_

> `Total_efficiency`	16.18	PEPPRO	_RES_

> `Read_depth`	3.41	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_1 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/tmpc6_9h7v5 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam` (14144,14149,14150)
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
[bam_sort_core] merging from 13 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:44:25. Running peak memory: 6.727GB.  
  PID: 14144;	Command: bowtie2;	Return code: 0;	Memory used: 3.764GB  
  PID: 14149;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 14150;	Command: samtools;	Return code: 0;	Memory used: 0.895GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam` (18832)
<pre>
</pre>
Command completed. Elapsed time: 0:01:29. Running peak memory: 6.727GB.  
  PID: 18832;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 09:44:09) elapsed: 6537.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R2.fq` (18927)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 6.727GB.  
  PID: 18927;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/prealignments/H9_treated_PRO-seq_1_human_rDNA_unmap_R1.fq` (18962)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 6.727GB.  
  PID: 18962;	Command: pigz;	Return code: 0;	Memory used: 0.011GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam` (19015)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 6.727GB.  
  PID: 19015;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	655778	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (19245)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.727GB.  
  PID: 19245;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/chr_sizes.bed` (19284,19285,19286,19287)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 19285;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 19287;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 19284;	Command: samtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 19286;	Command: awk;	Return code: 0;	Memory used: 0.0GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_noMT.bam` (19289)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 6.727GB.  
  PID: 19289;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (19347)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 19347;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam` (19348)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 6.727GB.  
  PID: 19348;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


### Split BAM file (06-15 09:47:05) elapsed: 176.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam` (19391,19392)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:35. Running peak memory: 6.727GB.  
  PID: 19391;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 19392;	Command: samtools;	Return code: 0;	Memory used: 5.36GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE2.bam` (19658,19659)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:02. Running peak memory: 6.727GB.  
  PID: 19658;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 19659;	Command: samtools;	Return code: 0;	Memory used: 4.161GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam` (20174)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 6.727GB.  
  PID: 20174;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam` (20231,20232)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:29. Running peak memory: 6.727GB.  
  PID: 20231;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20232;	Command: samtools;	Return code: 0;	Memory used: 5.071GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE2.bam` (20669,20670)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:01:55. Running peak memory: 6.727GB.  
  PID: 20669;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20670;	Command: samtools;	Return code: 0;	Memory used: 3.933GB


### Calculate library complexity (06-15 09:57:38) elapsed: 633.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam` (20851)
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
Command completed. Elapsed time: 0:01:41. Running peak memory: 6.727GB.  
  PID: 20851;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam` (20976)
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
............_.........................................................._......._.......................
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:53. Running peak memory: 6.727GB.  
  PID: 20976;	Command: preseq;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_counts.txt` (21425)
<pre>
psutil.NoSuchProcess process no longer exists (pid=21427)
Warning: couldn't add memory use for process: 21425
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.727GB.  
  PID: 21425;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_plot` (21453)
<pre>
Processing H9_treated_PRO-seq_1
INFO: Found real counts for H9_treated_PRO-seq_1 - Total (M): 19.73867 Unique (M): 18.664481

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.727GB.  
  PID: 21453;	Command: Rscript;	Return code: 0;	Memory used: 0.285GB

> `Library complexity`	QC_hg38/H9_treated_PRO-seq_1_preseq_plot.pdf	Library complexity	QC_hg38/H9_treated_PRO-seq_1_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.8094	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 10:01:45) elapsed: 247.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam` (21473)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.727GB.  
  PID: 21473;	Command: samtools;	Return code: 0;	Memory used: 0.016GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv` (21511)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/tmp_H9_treated_PRO-seq_1_PE1_bybpenda'
Processing with 12 cores...
Discarding 94 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000214v1']
Keeping 101 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 6.727GB.  
  PID: 21511;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 1.306GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	9869335.0	PEPPRO	_RES_

> `PBC2`	9869335.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_unmap.bam` (21562)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 6.727GB.  
  PID: 21562;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_temp.bam`

> `Unmapped_reads`	7715519	PEPPRO	_RES_

### Split BAM by strand (06-15 10:02:39) elapsed: 54.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam` (21605)
<pre>
</pre>
Command completed. Elapsed time: 0:01:05. Running peak memory: 6.727GB.  
  PID: 21605;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam` (21664)
<pre>
</pre>
Command completed. Elapsed time: 0:01:04. Running peak memory: 6.727GB.  
  PID: 21664;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 10:04:49) elapsed: 130.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (21719)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 21719;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt` (21721)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.727GB.  
  PID: 21721;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.981GB


> `TSS_coding_score`	53.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt` (21753)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.727GB.  
  PID: 21753;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.978GB


> `TSS_non-coding_score`	17.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt` (21880)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 6.727GB.  
  PID: 21880;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/H9_treated_PRO-seq_1_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_treated_PRO-seq_1_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt` (22033,22034,22035,22036)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22033;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 22035;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 22034;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 22036;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt` (22038)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22038;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 10:05:08) elapsed: 19.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed` (22040,22041)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.727GB.  
  PID: 22040;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 22041;	Command: bedtools;	Return code: 0;	Memory used: 0.095GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed` (22045,22046)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22045;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 22046;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSS_density.bed` (22048,22049,22050,22051)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 6.727GB.  
  PID: 22049;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 22051;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 22048;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 22050;	Command: sort;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_gene_body_density.bed` (22096,22097,22098)
<pre>
</pre>
Command completed. Elapsed time: 0:00:29. Running peak memory: 6.727GB.  
  PID: 22097;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 22096;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB  
  PID: 22098;	Command: sort;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/tmpunwmqyci` (22154,22155,22156)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22154;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 22156;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 22155;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/tmpunwmqyci | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0193886) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/tmpunwmqyci > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed` (22162)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22162;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	62.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed` (22167)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.727GB.  
  PID: 22167;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Pause index`	QC_hg38/H9_treated_PRO-seq_1_pause_index.pdf	Pause index	QC_hg38/H9_treated_PRO-seq_1_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_pause_index.bed` (22190)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22190;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 10:06:11) elapsed: 62.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam`
18009111.5 7282779

> `Plus_FRiP`	0.4	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam`
18009111.5 6888359

> `Minus_FRiP`	0.38	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_gene_sort.bed` (22316,22317)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 22316;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 22317;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_gene_coverage.bed` (22320)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 6.727GB.  
  PID: 22320;	Command: bedtools;	Return code: 0;	Memory used: 0.06GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed.gz` (22382)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22382;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed` (22383)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22383;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 10:07:19) elapsed: 68.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/raw/hg38_annotations.bed` (22392)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 22392;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed` (22394,22395,22396,22397)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 22394;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 22395;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 22397;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 22396;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed` (22400)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.727GB.  
  PID: 22400;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_minus_coverage.bed` (22413)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.727GB.  
  PID: 22413;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed` (22425,22426,22427,22428)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22425;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 22426;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 22428;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 22427;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed` (22430)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.727GB.  
  PID: 22430;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_minus_coverage.bed` (22445)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.727GB.  
  PID: 22445;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region"` (22457)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22457;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed` (22458,22459,22460,22461)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 22458;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 22460;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 22459;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 22461;	Command: bedtools;	Return code: 0;	Memory used: 0.051GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed` (22464)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.727GB.  
  PID: 22464;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_minus_coverage.bed` (22477)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.727GB.  
  PID: 22477;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR"` (22491)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22491;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed` (22492,22493,22494,22495)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 22492;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 22493;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 22495;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 22494;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed` (22497)
<pre>
</pre>
Command completed. Elapsed time: 0:00:13. Running peak memory: 6.727GB.  
  PID: 22497;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_minus_coverage.bed` (22510)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.727GB.  
  PID: 22510;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR"` (22522)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 22522;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed` (22523,22524,22525,22526)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 22523;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 22524;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 22526;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB  
  PID: 22525;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed` (22528)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.727GB.  
  PID: 22528;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_minus_coverage.bed` (22595)
<pre>
</pre>
Command completed. Elapsed time: 0:00:12. Running peak memory: 6.727GB.  
  PID: 22595;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed` (22642,22643,22644,22645)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 6.727GB.  
  PID: 22642;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 22643;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 22645;	Command: bedtools;	Return code: 0;	Memory used: 0.172GB  
  PID: 22644;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed` (22649)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.727GB.  
  PID: 22649;	Command: bedtools;	Return code: 0;	Memory used: 0.027GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_minus_coverage.bed` (22718)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.727GB.  
  PID: 22718;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed` (22770,22771,22772,22773)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 22770;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 22772;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 22771;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 22773;	Command: bedtools;	Return code: 0;	Memory used: 0.082GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed` (22776)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 6.727GB.  
  PID: 22776;	Command: bedtools;	Return code: 0;	Memory used: 0.076GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_minus_coverage.bed` (23057)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 6.727GB.  
  PID: 23057;	Command: bedtools;	Return code: 0;	Memory used: 0.038GB


### Plot cFRiF/FRiF (06-15 10:10:30) elapsed: 191.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_1 -z 3099922541 -n 10005525 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed` (23090)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 6.727GB.  
  PID: 23090;	Command: Rscript;	Return code: 0;	Memory used: 0.456GB

> `cFRiF`	QC_hg38/H9_treated_PRO-seq_1_cFRiF.pdf	cFRiF	QC_hg38/H9_treated_PRO-seq_1_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_1 -z 3099922541 -n 10005525 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_Intron_plus_coverage.bed` (23220)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.727GB.  
  PID: 23220;	Command: Rscript;	Return code: 0;	Memory used: 0.497GB

> `FRiF`	QC_hg38/H9_treated_PRO-seq_1_FRiF.pdf	FRiF	QC_hg38/H9_treated_PRO-seq_1_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 10:11:31) elapsed: 61.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_exons_sort.bed` (23262,23263)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.727GB.  
  PID: 23263;	Command: bedtools;	Return code: 0;	Memory used: 0.096GB  
  PID: 23262;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_introns_sort.bed` (23269,23270,23271)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.727GB.  
  PID: 23269;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 23271;	Command: bedtools;	Return code: 0;	Memory used: 0.003GB  
  PID: 23270;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_coverage.bed` (23283)
<pre>
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 6.727GB.  
  PID: 23283;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_coverage.bed` (23323)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 6.727GB.  
  PID: 23323;	Command: bedtools;	Return code: 0;	Memory used: 0.078GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.0091115)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_rpkm.bed` (23375,23376,23377)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 23375;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 23377;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 23376;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/18.0091115)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_rpkm.bed` (23380,23381,23382)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.727GB.  
  PID: 23380;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 23382;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 23381;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed` (23384,23385,23386)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 23384;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 23386;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 23385;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.19	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed --annotate` (23392)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 6.727GB.  
  PID: 23392;	Command: Rscript;	Return code: 0;	Memory used: 0.319GB

> `mRNA contamination`	QC_hg38/H9_treated_PRO-seq_1_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_treated_PRO-seq_1_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/QC_hg38/H9_treated_PRO-seq_1_exon_intron_ratios.bed` (23420)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.727GB.  
  PID: 23420;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (06-15 10:12:39) elapsed: 69.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam` (23428)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.727GB.  
  PID: 23428;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 18009111.5` (23437)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_plus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_1_plus_cuttrace_mzrp_pq2'
Processing with 4 cores...
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270720v1_random', 'chr22_KI270735v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_GL000214v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270752v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:18. Running peak memory: 6.727GB.  
  PID: 23437;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.785GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam` (25650)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 6.727GB.  
  PID: 25650;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 18009111.5` (25735)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/aligned_hg38/H9_treated_PRO-seq_1_minus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_1_minus_cuttrace_evfj1ot8'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 106 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr17_KI270730v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1']
Keeping 89 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 89 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_exact_body_0-mer.bw'
Merging 89 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_1/signal_hg38/H9_treated_PRO-seq_1_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:17. Running peak memory: 6.727GB.  
  PID: 25735;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.731GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:10:12
*  Total elapsed time (all runs):  5:46:27
*         Peak memory (this run):  6.7273 GB
*        Pipeline completed time: 2020-06-15 10:27:30
