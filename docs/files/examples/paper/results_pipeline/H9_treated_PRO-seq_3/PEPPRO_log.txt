### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name H9_treated_PRO-seq_3 --genome hg38 --input /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep3_PE1.fastq.gz --single-or-paired PAIRED -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 12 -M 12000 --input2 /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep3_PE2.fastq.gz --protocol PRO --umi-len 8 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-28c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/
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
*              `input`:  `['/project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep3_PE1.fastq.gz']`
*             `input2`:  `['/project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep3_PE2.fastq.gz']`
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
*        `sample_name`:  `H9_treated_PRO-seq_3`
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

Local input file: /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep3_PE1.fastq.gz
Local input file: /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep3_PE2.fastq.gz

> `File_mb`	2905.11	PEPPRO	_RES_

> `Read_type`	PAIRED	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-15 07:17:19) elapsed: 1.0 _TIME_

Number of input file sets: 2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/H9_treated_PRO-seq_3_R1.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep3_PE1.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/H9_treated_PRO-seq_3_R1.fastq.gz` (138710)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 138710;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/H9_treated_PRO-seq_3_R1.fastq.gz'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/H9_treated_PRO-seq_3_R2.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/H9_200nM_romidepsin_rep3_PE2.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/H9_treated_PRO-seq_3_R2.fastq.gz` (138713)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 138713;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/H9_treated_PRO-seq_3_R2.fastq.gz'
Found .fastq.gz file
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1.fastq`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2.fastq`  

> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/H9_treated_PRO-seq_3_R1.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1.fastq` (138720)
<pre>
</pre>
Command completed. Elapsed time: 0:02:52. Running peak memory: 0.002GB.  
  PID: 138720;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/H9_treated_PRO-seq_3_R2.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2.fastq` (54553)
<pre>
</pre>
Command completed. Elapsed time: 0:01:56. Running peak memory: 0.002GB.  
  PID: 54553;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


> `Raw_reads`	120377388	PEPPRO	_RES_

> `Fastq_reads`	120377388	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/H9_treated_PRO-seq_3_R1.fastq.gz', '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/H9_treated_PRO-seq_3_R2.fastq.gz']

### FASTQ processing:  (06-15 07:25:41) elapsed: 502.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3_R1_cutadapt.txt` (3572)
<pre>
</pre>
Command completed. Elapsed time: 0:02:51. Running peak memory: 4.399GB.  
  PID: 3572;	Command: cutadapt;	Return code: 0;	Memory used: 4.399GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq` (3798,3799)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 4.399GB.  
  PID: 3798;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 3799;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	32320459	PEPPRO	_RES_

> `Trim_loss_rate`	73.15	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq` (3889)
<pre>
Started analysis of H9_treated_PRO-seq_3_R1_processed.fastq
Approx 5% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 10% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 15% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 20% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 25% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 30% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 35% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 40% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 45% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 50% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 55% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 60% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 65% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 70% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 75% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 80% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 85% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 90% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 95% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Analysis complete for H9_treated_PRO-seq_3_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:03. Running peak memory: 4.399GB.  
  PID: 3889;	Command: fastqc;	Return code: 0;	Memory used: 0.178GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_3_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_trimmed.fastq`  

> `seqkit rmdup --threads 12 --by-seq --ignore-case -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_dedup.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_noadap.fastq` (4233)
<pre>
[INFO][0m 4932446 duplicated records removed
</pre>
Command completed. Elapsed time: 0:01:20. Running peak memory: 4.399GB.  
  PID: 4233;	Command: seqkit;	Return code: 0;	Memory used: 2.046GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_dedup.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_trimmed.fastq` (4357,4358)
<pre>
</pre>
Command completed. Elapsed time: 0:00:46. Running peak memory: 4.399GB.  
  PID: 4357;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 4358;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	46004533.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	26667722.0	PEPPRO	_RES_

> `Duplicate_reads`	4932446.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	44.3069	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/processed_R1.flag` (4567)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 4.399GB.  
  PID: 4567;	Command: touch;	Return code: 0;	Memory used: 0.002GB


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed.fastq`  

> `(cutadapt -j 12 -m 10 -O 1 -a GATCGTCGGACTGTAGAACTCTGAAC /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3_R2_cutadapt.txt` (4570)
<pre>
</pre>
Command completed. Elapsed time: 0:02:00. Running peak memory: 4.399GB.  
  PID: 4570;	Command: cutadapt;	Return code: 0;	Memory used: 3.415GB


> `seqtk trimfq -b 8 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_noadap.fastq | seqtk seq -L 10 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed.fastq` (4947,4948)
<pre>
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 4.399GB.  
  PID: 4947;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 4948;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB

Evaluating read trimming

> `Trimmed_reads`	64640918	PEPPRO	_RES_

> `Trim_loss_rate`	46.3	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq` (5051)
<pre>
Started analysis of H9_treated_PRO-seq_3_R1_processed.fastq
Approx 5% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 10% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 15% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 20% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 25% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 30% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 35% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 40% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 45% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 50% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 55% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 60% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 65% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 70% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 75% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 80% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 85% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 90% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Approx 95% complete for H9_treated_PRO-seq_3_R1_processed.fastq
Analysis complete for H9_treated_PRO-seq_3_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:01:02. Running peak memory: 4.399GB.  
  PID: 5051;	Command: fastqc;	Return code: 0;	Memory used: 0.179GB

> `FastQC report r1`	fastqc/H9_treated_PRO-seq_3_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed.fastq` (5140)
<pre>
Started analysis of H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 5% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 10% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 15% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 20% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 25% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 30% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 35% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 40% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 45% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 50% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 55% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 60% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 65% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 70% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 75% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 80% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 85% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 90% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Approx 95% complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
Analysis complete for H9_treated_PRO-seq_3_R2_trimmed.fastq
</pre>
Command completed. Elapsed time: 0:01:19. Running peak memory: 4.399GB.  
  PID: 5140;	Command: fastqc;	Return code: 0;	Memory used: 0.165GB

> `FastQC report r2`	fastqc/H9_treated_PRO-seq_3_R2_trimmed_fastqc.html	FastQC report r2	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3.hist`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3.histogram`  

> `fastq_pair -t 108339649 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_noadap.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_noadap.fastq` (5227)
<pre>
Left paired: 33366507		Right paired: 33366507
Left single: 154465		Right single: 1916320
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_noadap.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_noadap.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_noadap.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_noadap.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:04:36. Running peak memory: 7.355GB.  
  PID: 5227;	Command: fastq_pair;	Return code: 0;	Memory used: 7.355GB


> `flash -q -t 12 --compress-prog=pigz --suffix=gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_noadap.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_noadap.fastq.paired.fq -o H9_treated_PRO-seq_3 -d /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt` (6033)
<pre>
</pre>
Command completed. Elapsed time: 0:01:22. Running peak memory: 7.355GB.  
  PID: 6033;	Command: flash;	Return code: 0;	Memory used: 0.091GB


### Plot adapter insertion distribution (06-15 07:45:31) elapsed: 1189.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R adapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3.hist -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt -u 8` (6460)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 7.355GB.  
  PID: 6460;	Command: Rscript;	Return code: 0;	Memory used: 0.203GB

> `Adapter insertion distribution`	cutadapt/H9_treated_PRO-seq_3_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/H9_treated_PRO-seq_3_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3.hist | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print len-8}'`

> `Peak_adapter_insertion_size`	20	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-15 07:45:37) elapsed: 6.0 _TIME_


> `awk '{ if (($1-8) == 10) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 20) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 30) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3.hist`

> `awk '{ if (($1-8) == 40) {status = 1}} END {if (status) {print status} else {print 0}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3.hist`

> `awk '(($1-8 ) <= 20 && ($1-8 ) >= 10){degradedSum += $2}; (($1-8 ) >= 30 && ($1-8) <= 40){intactSum += $2}  END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/cutadapt/H9_treated_PRO-seq_3.hist`

> `Degradation_ratio`	0.9116	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed_dups.fastq`  

> `cp /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed_dups.fastq` (6490)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 7.355GB.  
  PID: 6490;	Command: cp;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/processed_R2.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/processed_R2.flag` (6516)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 6516;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/repaired.flag`  

> `fastq_pair -t 108339649 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed.fastq` (6517)
<pre>
Left paired: 32136774		Right paired: 32136774
Left single: 183685		Right single: 1977423
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:03:57. Running peak memory: 7.355GB.  
  PID: 6517;	Command: fastq_pair;	Return code: 0;	Memory used: 7.13GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq` (6742)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 6742;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed.fastq` (6743)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 6743;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/repaired.flag` (6745)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 6745;	Command: touch;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/dups_repaired.flag`  

> `fastq_pair -t 108339649 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed_dups.fastq` (6746)
<pre>
Left paired: 27441885		Right paired: 27441885
Left single: 119928		Right single: 6672312
Writing the paired reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_trimmed.fastq.paired.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed_dups.fastq.paired.fq.
Writing the single reads to /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_trimmed.fastq.single.fq and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed_dups.fastq.single.fq
</pre>
Command completed. Elapsed time: 0:02:26. Running peak memory: 7.355GB.  
  PID: 6746;	Command: fastq_pair;	Return code: 0;	Memory used: 5.969GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_trimmed.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_trimmed.fastq` (7165)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 7165;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed_dups.fastq.paired.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed_dups.fastq` (7166)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 7166;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/dups_repaired.flag` (7168)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 7168;	Command: touch;	Return code: 0;	Memory used: 0.0GB


### Prealignments (06-15 07:52:24) elapsed: 407.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-15 07:52:24) elapsed: 0.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/human_rDNA_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/human_rDNA_bt2` (7169)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 7169;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/human_rDNA_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_R2.fq` (7170)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_3 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/human_rDNA_bt2 2>&1 > /dev/null)`
not gzipping output
File not added to cleanup: prealignments/human_rDNA_bt2
Missing stat 'Aligned_reads_human_rDNA'
32136774 reads; of these:
  32136774 (100.00%) were unpaired; of these:
    28788351 (89.58%) aligned 0 times
    3348423 (10.42%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
10.42% overall alignment rate

> `Aligned_reads_human_rDNA`	6696846.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	10.36	PEPPRO	_RES_

### Map to human_rDNA (06-15 07:55:57) elapsed: 214.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/human_rDNA_dups_bt2`  

> `mkfifo /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/human_rDNA_dups_bt2` (7602)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 7602;	Command: mkfifo;	Return code: 0;	Memory used: 0.0GB

File not added to cleanup: prealignments/human_rDNA_dups_bt2
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `perl /scratch/jps3dp/tools/databio//peppro/tools/filter_paired_fq.pl /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/human_rDNA_dups_bt2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_trimmed.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R2_trimmed_dups.fastq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_dups_R1.fq /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_dups_R2.fq` (7603)
<pre>
</pre>

> `(bowtie2 -p 12 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id H9_treated_PRO-seq_3 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/fastq/H9_treated_PRO-seq_3_R1_trimmed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/human_rDNA_dups_bt2 2>&1 > /dev/null)`
not gzipping output
3348423 reads skipped
0 reads lost
File not added to cleanup: prealignments/human_rDNA_dups_bt2

### Map to genome (06-15 07:59:20) elapsed: 203.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_fail_qc_dups.bam
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_3 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/tmpm93otle_ -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp.bam` (7809,7810,7811)
<pre>
2663663 reads skipped
0 reads lost
28788351 reads; of these:
  28788351 (100.00%) were paired; of these:
    12578156 (43.69%) aligned concordantly 0 times
    13526777 (46.99%) aligned concordantly exactly 1 time
    2683418 (9.32%) aligned concordantly >1 times
    ----
    12578156 pairs aligned concordantly 0 times; of these:
      3397179 (27.01%) aligned discordantly 1 time
    ----
    9180977 pairs aligned 0 times concordantly or discordantly; of these:
      18361954 mates make up the pairs; of these:
        7400125 (40.30%) aligned 0 times
        4038486 (21.99%) aligned exactly 1 time
        6923343 (37.70%) aligned >1 times
87.15% overall alignment rate
[bam_sort_core] merging from 16 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:48:46. Running peak memory: 7.355GB.  
  PID: 7810;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 7809;	Command: bowtie2;	Return code: 0;	Memory used: 3.772GB  
  PID: 7811;	Command: samtools;	Return code: 0;	Memory used: 0.895GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort.bam` (13121)
<pre>
</pre>
Command completed. Elapsed time: 0:01:35. Running peak memory: 7.355GB.  
  PID: 13121;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	50176577	PEPPRO	_RES_

> `QC_filtered_reads`	29080016	PEPPRO	_RES_

> `Aligned_reads`	21096561.0	PEPPRO	_RES_

> `Alignment_rate`	32.64	PEPPRO	_RES_

> `Total_efficiency`	17.53	PEPPRO	_RES_

> `Read_depth`	3.98	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort_dups.bam`  

> `bowtie2 -p 12 --very-sensitive -X 2000 --rg-id H9_treated_PRO-seq_3 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 --rf -1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_dups_R1.fq -2 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_dups_R2.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/tmpm93otle_ -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp_dups.bam` (14710,14711,14712)
<pre>
24778222 reads; of these:
  24778222 (100.00%) were paired; of these:
    10669067 (43.06%) aligned concordantly 0 times
    11790699 (47.58%) aligned concordantly exactly 1 time
    2318456 (9.36%) aligned concordantly >1 times
    ----
    10669067 pairs aligned concordantly 0 times; of these:
      2974174 (27.88%) aligned discordantly 1 time
    ----
    7694893 pairs aligned 0 times concordantly or discordantly; of these:
      15389786 mates make up the pairs; of these:
        6244577 (40.58%) aligned 0 times
        3533931 (22.96%) aligned exactly 1 time
        5611278 (36.46%) aligned >1 times
87.40% overall alignment rate
[bam_sort_core] merging from 14 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:42:51. Running peak memory: 7.355GB.  
  PID: 14710;	Command: bowtie2;	Return code: 0;	Memory used: 3.758GB  
  PID: 14711;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 14712;	Command: samtools;	Return code: 0;	Memory used: 0.895GB


> `samtools view -q 10 -b -@ 12 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_fail_qc_dups.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp_dups.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort_dups.bam` (18989)
<pre>
</pre>
Command completed. Elapsed time: 0:01:26. Running peak memory: 7.355GB.  
  PID: 18989;	Command: samtools;	Return code: 0;	Memory used: 0.018GB


### Compress all unmapped read files (06-15 09:46:12) elapsed: 6412.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_R2.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_R2.fq` (19314)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 7.355GB.  
  PID: 19314;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_R1.fq.gz`  

> `pigz -f -p 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/prealignments/H9_treated_PRO-seq_3_human_rDNA_unmap_R1.fq` (19353)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 7.355GB.  
  PID: 19353;	Command: pigz;	Return code: 0;	Memory used: 0.009GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp.bam` (19407)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 7.355GB.  
  PID: 19407;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	667197	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort.bam` (19473)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 7.355GB.  
  PID: 19473;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/chr_sizes.bed` (19501,19502,19503,19504)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 19503;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 19501;	Command: samtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 19504;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 19502;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/chr_sizes.bed -b -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_noMT.bam` (19506)
<pre>
</pre>
Command completed. Elapsed time: 0:00:39. Running peak memory: 7.355GB.  
  PID: 19506;	Command: samtools;	Return code: 0;	Memory used: 0.019GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort.bam` (19559)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 19559;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort.bam` (19560)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 7.355GB.  
  PID: 19560;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


### Split BAM file (06-15 09:49:29) elapsed: 197.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam` (19633,19634)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:53. Running peak memory: 7.355GB.  
  PID: 19633;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 19634;	Command: samtools;	Return code: 0;	Memory used: 6.208GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE2.bam` (20143,20144)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:23. Running peak memory: 7.355GB.  
  PID: 20144;	Command: samtools;	Return code: 0;	Memory used: 4.661GB  
  PID: 20143;	Command: samtools;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	38	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp_dups.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp_dups.bam` (20692)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 7.355GB.  
  PID: 20692;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp_dups.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort_dups.bam`  
No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort_dups.bam.bai
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_dups_PE1.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_dups_PE2.bam`  

> `samtools view -b -f 64 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_dups_PE1.bam` (20726,20727)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:36. Running peak memory: 7.355GB.  
  PID: 20726;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20727;	Command: samtools;	Return code: 0;	Memory used: 5.495GB


> `samtools view -b -f 128 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_sort_dups.bam | samtools sort - -@ 12 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_dups_PE2.bam` (20945,20946)
<pre>
[bam_sort_core] merging from 0 files and 12 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:02:07. Running peak memory: 7.355GB.  
  PID: 20945;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 20946;	Command: samtools;	Return code: 0;	Memory used: 4.391GB


### Calculate library complexity (06-15 10:01:09) elapsed: 700.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_out.txt`  

> `preseq c_curve -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_out.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_dups_PE1.bam` (21422)
<pre>
BAM_INPUT
TOTAL READS     = 20168929
COUNTS_SUM      = 20168929
DISTINCT READS  = 1.4613e+07
DISTINCT COUNTS = 320
MAX COUNT       = 14961
COUNTS OF 1     = 1.25301e+07
OBSERVED COUNTS (14962)
1	12530089
2	1265230
3	339731
4	151204
5	85421
6	54520
7	37501
8	26746
9	19937
10	15846
11	12083
12	9949
13	8065
14	6771
15	5607
16	4757
17	4010
18	3469
19	2981
20	2653
21	2307
22	1993
23	1854
24	1591
25	1387
26	1247
27	1131
28	985
29	895
30	810
31	746
32	697
33	666
34	606
35	568
36	521
37	480
38	454
39	380
40	365
41	350
42	333
43	308
44	273
45	240
46	229
47	211
48	234
49	223
50	197
51	159
52	172
53	144
54	136
55	128
56	137
57	129
58	129
59	122
60	112
61	86
62	95
63	83
64	99
65	73
66	77
67	71
68	63
69	73
70	69
71	60
72	57
73	57
74	50
75	66
76	49
77	41
78	48
79	42
80	39
81	31
82	43
83	43
84	45
85	34
86	32
87	36
88	34
89	21
90	24
91	31
92	34
93	30
94	27
95	31
96	19
97	24
98	17
99	24
100	19
101	19
102	18
103	12
104	25
105	15
106	13
107	15
108	19
109	17
110	16
111	14
112	13
113	15
114	14
115	9
116	17
117	10
118	16
119	18
120	14
121	10
122	12
123	18
124	6
125	9
126	8
127	10
128	7
129	5
130	6
131	8
132	11
133	6
134	10
135	9
136	5
137	6
138	8
139	2
140	8
141	7
142	4
143	3
144	7
145	7
146	6
147	10
148	9
149	5
150	8
151	4
152	7
153	5
154	1
155	8
156	3
157	3
158	4
159	7
160	2
161	7
162	4
163	5
164	5
165	2
166	2
167	3
168	4
169	4
170	1
171	4
172	4
173	2
174	5
175	2
176	1
177	3
178	2
179	1
180	4
181	1
182	1
183	5
184	2
185	3
186	2
187	3
188	4
189	2
190	1
191	4
192	2
193	2
196	3
197	4
199	2
200	1
201	2
202	1
203	1
204	1
205	1
206	1
207	1
209	1
210	2
211	1
212	1
213	3
214	3
215	2
216	3
217	1
218	1
219	1
222	2
223	1
224	1
225	1
226	1
229	1
230	1
231	4
232	3
233	1
237	2
238	2
239	2
241	1
242	1
244	3
246	1
247	1
251	2
254	2
255	3
256	3
259	2
260	1
262	1
264	2
265	1
268	2
269	1
270	1
271	1
273	1
276	1
278	3
280	1
284	1
289	1
291	3
295	1
302	1
303	1
305	2
306	1
308	1
309	1
311	1
314	2
315	2
319	1
325	2
328	1
329	1
330	1
340	2
342	1
343	1
346	1
347	1
355	1
377	1
378	1
379	1
389	1
390	1
396	1
400	1
404	1
405	1
406	1
410	1
416	1
421	1
426	1
427	1
432	1
441	1
453	1
454	1
459	1
461	1
504	1
512	1
514	1
551	1
598	1
606	1
646	1
686	1
693	1
727	1
736	1
792	1
858	1
903	1
951	1
1114	1
1290	1
1291	1
1482	1
1554	1
3863	1
4886	1
5882	1
12486	1
14961	1

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
sample size: 20000000
</pre>
Command completed. Elapsed time: 0:01:52. Running peak memory: 7.355GB.  
  PID: 21422;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_yield.txt`  

> `preseq lc_extrap -v -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_yield.txt -B /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_dups_PE1.bam` (21625)
<pre>
BAM_INPUT
TOTAL READS     = 20168929
DISTINCT READS  = 1.4613e+07
DISTINCT COUNTS = 320
MAX COUNT       = 14961
COUNTS OF 1     = 1.25301e+07
MAX TERMS       = 100
OBSERVED COUNTS (14962)
1	12530089
2	1265230
3	339731
4	151204
5	85421
6	54520
7	37501
8	26746
9	19937
10	15846
11	12083
12	9949
13	8065
14	6771
15	5607
16	4757
17	4010
18	3469
19	2981
20	2653
21	2307
22	1993
23	1854
24	1591
25	1387
26	1247
27	1131
28	985
29	895
30	810
31	746
32	697
33	666
34	606
35	568
36	521
37	480
38	454
39	380
40	365
41	350
42	333
43	308
44	273
45	240
46	229
47	211
48	234
49	223
50	197
51	159
52	172
53	144
54	136
55	128
56	137
57	129
58	129
59	122
60	112
61	86
62	95
63	83
64	99
65	73
66	77
67	71
68	63
69	73
70	69
71	60
72	57
73	57
74	50
75	66
76	49
77	41
78	48
79	42
80	39
81	31
82	43
83	43
84	45
85	34
86	32
87	36
88	34
89	21
90	24
91	31
92	34
93	30
94	27
95	31
96	19
97	24
98	17
99	24
100	19
101	19
102	18
103	12
104	25
105	15
106	13
107	15
108	19
109	17
110	16
111	14
112	13
113	15
114	14
115	9
116	17
117	10
118	16
119	18
120	14
121	10
122	12
123	18
124	6
125	9
126	8
127	10
128	7
129	5
130	6
131	8
132	11
133	6
134	10
135	9
136	5
137	6
138	8
139	2
140	8
141	7
142	4
143	3
144	7
145	7
146	6
147	10
148	9
149	5
150	8
151	4
152	7
153	5
154	1
155	8
156	3
157	3
158	4
159	7
160	2
161	7
162	4
163	5
164	5
165	2
166	2
167	3
168	4
169	4
170	1
171	4
172	4
173	2
174	5
175	2
176	1
177	3
178	2
179	1
180	4
181	1
182	1
183	5
184	2
185	3
186	2
187	3
188	4
189	2
190	1
191	4
192	2
193	2
196	3
197	4
199	2
200	1
201	2
202	1
203	1
204	1
205	1
206	1
207	1
209	1
210	2
211	1
212	1
213	3
214	3
215	2
216	3
217	1
218	1
219	1
222	2
223	1
224	1
225	1
226	1
229	1
230	1
231	4
232	3
233	1
237	2
238	2
239	2
241	1
242	1
244	3
246	1
247	1
251	2
254	2
255	3
256	3
259	2
260	1
262	1
264	2
265	1
268	2
269	1
270	1
271	1
273	1
276	1
278	3
280	1
284	1
289	1
291	3
295	1
302	1
303	1
305	2
306	1
308	1
309	1
311	1
314	2
315	2
319	1
325	2
328	1
329	1
330	1
340	2
342	1
343	1
346	1
347	1
355	1
377	1
378	1
379	1
389	1
390	1
396	1
400	1
404	1
405	1
406	1
410	1
416	1
421	1
426	1
427	1
432	1
441	1
453	1
454	1
459	1
461	1
504	1
512	1
514	1
551	1
598	1
606	1
646	1
686	1
693	1
727	1
736	1
792	1
858	1
903	1
951	1
1114	1
1290	1
1291	1
1482	1
1554	1
3863	1
4886	1
5882	1
12486	1
14961	1

[ESTIMATING YIELD CURVE]
[BOOTSTRAPPING HISTOGRAM]
..............................................................................................._.....
[COMPUTING CONFIDENCE INTERVALS]
[WRITING OUTPUT]
</pre>
Command completed. Elapsed time: 0:01:57. Running peak memory: 7.355GB.  
  PID: 21625;	Command: preseq;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_counts.txt`  

> `echo '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_yield.txt '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_dups_PE1.bam)' '$(samtools view -c -F 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_counts.txt` (21786)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 7.355GB.  
  PID: 21786;	Command: echo;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_plot.pdf`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_plot.png`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R preseq -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_yield.txt -r /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_counts.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_plot` (22074)
<pre>
Processing H9_treated_PRO-seq_3
INFO: Found real counts for H9_treated_PRO-seq_3 - Total (M): 22.825749 Unique (M): 20.168929

Library complexity plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 7.355GB.  
  PID: 22074;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `Library complexity`	QC_hg38/H9_treated_PRO-seq_3_preseq_plot.pdf	Library complexity	QC_hg38/H9_treated_PRO-seq_3_preseq_plot.png	PEPPRO	_OBJ_
Missing stat 'Frac_exp_unique_at_10M'

> `grep -w '10000000' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_preseq_yield.txt | awk '{print $2}'`

> `Frac_exp_unique_at_10M`	0.7916	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-15 10:05:33) elapsed: 264.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam` (22092)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 7.355GB.  
  PID: 22092;	Command: samtools;	Return code: 0;	Memory used: 0.015GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam -c 12 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_bamQC.tsv` (22112)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/tmp_H9_treated_PRO-seq_3_PE1_89bj3hev'
Processing with 12 cores...
Discarding 98 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1']
Keeping 97 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:00:24. Running peak memory: 7.355GB.  
  PID: 22112;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 2.007GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_bamQC.tsv`

> `NRF`	1.0	PEPPRO	_RES_

> `PBC1`	11412874.5	PEPPRO	_RES_

> `PBC2`	11412874.5	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_unmap.bam`  

> `samtools view -b -@ 12 -f 12  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_unmap.bam` (22208)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 7.355GB.  
  PID: 22208;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools view -c -f 4 -@ 12 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_temp.bam`

> `Unmapped_reads`	7400125	PEPPRO	_RES_

### Split BAM by strand (06-15 10:06:33) elapsed: 60.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam` (22277)
<pre>
</pre>
Command completed. Elapsed time: 0:01:17. Running peak memory: 7.355GB.  
  PID: 22277;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam` (22436)
<pre>
</pre>
Command completed. Elapsed time: 0:01:15. Running peak memory: 7.355GB.  
  PID: 22436;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-15 10:09:05) elapsed: 152.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (22531)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 22531;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/plus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_plus_TssEnrichment.txt` (22532)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 7.355GB.  
  PID: 22532;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.702GB


> `TSS_coding_score`	60.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/minus_TSS.tsv -p ends -c 12 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_minus_TssEnrichment.txt` (22565)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 7.355GB.  
  PID: 22565;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.767GB


> `TSS_non-coding_score`	17.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_minus_TssEnrichment.txt` (22599)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 7.355GB.  
  PID: 22599;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/H9_treated_PRO-seq_3_TSSenrichment.pdf	TSS enrichment	QC_hg38/H9_treated_PRO-seq_3_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt` (22622,22623,22624,22625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 22622;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 22624;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 22623;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 22625;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_keep.txt` (22627)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 22627;	Command: cut;	Return code: 0;	Memory used: 0.0GB


### Calculate Pause Index (PI) (06-15 10:09:26) elapsed: 21.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed` (22629,22630)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 7.355GB.  
  PID: 22629;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 22630;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed` (22633,22634)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 22633;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 22634;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_TSS_density.bed` (22636,22637,22638,22639)
<pre>
</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 7.355GB.  
  PID: 22637;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 22639;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 22636;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB  
  PID: 22638;	Command: sort;	Return code: 0;	Memory used: 0.004GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_gene_body_density.bed` (22766,22767,22768)
<pre>
</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 7.355GB.  
  PID: 22768;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 22766;	Command: bedtools;	Return code: 0;	Memory used: 0.057GB  
  PID: 22767;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/tmpq0v5ymlg` (23104,23105,23106)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 23104;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 23106;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 23105;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/tmpq0v5ymlg | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.0223084) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/tmpq0v5ymlg > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_pause_index.bed` (23112)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 23112;	Command: awk;	Return code: 0;	Memory used: 0.002GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	63.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_pause_index.bed` (23117)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 7.355GB.  
  PID: 23117;	Command: Rscript;	Return code: 0;	Memory used: 0.317GB

> `Pause index`	QC_hg38/H9_treated_PRO-seq_3_pause_index.pdf	Pause index	QC_hg38/H9_treated_PRO-seq_3_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_pause_index.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_pause_index.bed` (23138)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 23138;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-15 10:10:40) elapsed: 74.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam`
21096561.0 8414821

> `Plus_FRiP`	0.4	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam`
21096561.0 7977228

> `Minus_FRiP`	0.38	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_gene_sort.bed` (23249,23250)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 23249;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 23250;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_gene_coverage.bed` (23252)
<pre>
</pre>
Command completed. Elapsed time: 0:00:35. Running peak memory: 7.355GB.  
  PID: 23252;	Command: bedtools;	Return code: 0;	Memory used: 0.06GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/hg38_annotations.bed.gz` (23297)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 23297;	Command: ln;	Return code: 0;	Memory used: 0.0GB


> `pigz -f -p 12 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/hg38_annotations.bed` (23298)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 23298;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-15 10:11:57) elapsed: 77.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/raw/hg38_annotations.bed` (23308)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 23308;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Enhancer_sort.bed` (23310,23311,23312,23313)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 23310;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23311;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 23313;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 23312;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Enhancer_plus_coverage.bed` (23316)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 7.355GB.  
  PID: 23316;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Enhancer_minus_coverage.bed` (23349)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 7.355GB.  
  PID: 23349;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_sort.bed` (23365,23366,23367,23368)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 23365;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23366;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 23368;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 23367;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_plus_coverage.bed` (23370)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 7.355GB.  
  PID: 23370;	Command: bedtools;	Return code: 0;	Memory used: 0.023GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_minus_coverage.bed` (23433)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 7.355GB.  
  PID: 23433;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_Flanking_Region"` (23508)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 23508;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed` (23509,23510,23511,23512)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 23509;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23511;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 23510;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 23512;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed` (23514)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 7.355GB.  
  PID: 23514;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_Flanking_Region_minus_coverage.bed` (23545)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 7.355GB.  
  PID: 23545;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/5_UTR"` (23572)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 23572;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/5_UTR_sort.bed` (23573,23574,23575,23576)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 23573;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23574;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 23576;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 23575;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_5_UTR_plus_coverage.bed` (23578)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 7.355GB.  
  PID: 23578;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_5_UTR_minus_coverage.bed` (23627)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 7.355GB.  
  PID: 23627;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/3_UTR"` (23669)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 23669;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/3_UTR_sort.bed` (23670,23671,23672,23673)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 23670;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23671;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 23673;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB  
  PID: 23672;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_3_UTR_plus_coverage.bed` (23675)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 7.355GB.  
  PID: 23675;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_3_UTR_minus_coverage.bed` (23701)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 7.355GB.  
  PID: 23701;	Command: bedtools;	Return code: 0;	Memory used: 0.011GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Exon_sort.bed` (23739,23740,23741,23742)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 7.355GB.  
  PID: 23739;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23740;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 23742;	Command: bedtools;	Return code: 0;	Memory used: 0.158GB  
  PID: 23741;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Exon_plus_coverage.bed` (23746)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 7.355GB.  
  PID: 23746;	Command: bedtools;	Return code: 0;	Memory used: 0.025GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Exon_minus_coverage.bed` (23782)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 7.355GB.  
  PID: 23782;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Intron_sort.bed` (23995,23996,23997,23998)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 23995;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 23997;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 23996;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 23998;	Command: bedtools;	Return code: 0;	Memory used: 0.077GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Intron_plus_coverage.bed` (24046)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 7.355GB.  
  PID: 24046;	Command: bedtools;	Return code: 0;	Memory used: 0.058GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Intron_minus_coverage.bed` (24598)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 7.355GB.  
  PID: 24598;	Command: bedtools;	Return code: 0;	Memory used: 0.037GB


### Plot cFRiF/FRiF (06-15 10:15:40) elapsed: 222.0 _TIME_


> `samtools view -@ 12 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_3 -z 3099922541 -n 11561442 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Intron_plus_coverage.bed` (24663)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 7.355GB.  
  PID: 24663;	Command: Rscript;	Return code: 0;	Memory used: 0.44GB

> `cFRiF`	QC_hg38/H9_treated_PRO-seq_3_cFRiF.pdf	cFRiF	QC_hg38/H9_treated_PRO-seq_3_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s H9_treated_PRO-seq_3 -z 3099922541 -n 11561442 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_Intron_plus_coverage.bed` (24726)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 7.355GB.  
  PID: 24726;	Command: Rscript;	Return code: 0;	Memory used: 0.448GB

> `FRiF`	QC_hg38/H9_treated_PRO-seq_3_FRiF.pdf	FRiF	QC_hg38/H9_treated_PRO-seq_3_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-15 10:16:43) elapsed: 64.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_exons_sort.bed` (24785,24786)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 7.355GB.  
  PID: 24786;	Command: bedtools;	Return code: 0;	Memory used: 0.089GB  
  PID: 24785;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_introns_sort.bed` (24794,24795,24796)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 7.355GB.  
  PID: 24794;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 24796;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 24795;	Command: bedtools;	Return code: 0;	Memory used: 0.033GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exons_coverage.bed` (24802)
<pre>
</pre>
Command completed. Elapsed time: 0:00:28. Running peak memory: 7.355GB.  
  PID: 24802;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_PE1.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_introns_coverage.bed` (24831)
<pre>
</pre>
Command completed. Elapsed time: 0:00:33. Running peak memory: 7.355GB.  
  PID: 24831;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/21.096561)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exons_rpkm.bed` (24968,24969,24970)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 24968;	Command: awk;	Return code: 0;	Memory used: 0.009GB  
  PID: 24970;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 24969;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/21.096561)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_introns_rpkm.bed` (24972,24973,24974)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 7.355GB.  
  PID: 24972;	Command: awk;	Return code: 0;	Memory used: 0.008GB  
  PID: 24974;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 24973;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exon_intron_ratios.bed` (24977,24978,24979)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 24977;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 24979;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 24978;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.1	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exon_intron_ratios.bed --annotate` (24985)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 7.355GB.  
  PID: 24985;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `mRNA contamination`	QC_hg38/H9_treated_PRO-seq_3_mRNA_contamination.pdf	mRNA contamination	QC_hg38/H9_treated_PRO-seq_3_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exon_intron_ratios.bed.gz`  

> `pigz -f -p 12 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/QC_hg38/H9_treated_PRO-seq_3_exon_intron_ratios.bed` (25006)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 7.355GB.  
  PID: 25006;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-15 10:18:01) elapsed: 78.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam` (25014)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 7.355GB.  
  PID: 25014;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_plus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 21096561.0` (25022)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_plus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_3_plus_cuttrace_mrjp2wis'
Processing with 4 cores...
stdin is empty of data
stdin is empty of data
Discarding 114 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr22_KI270735v1_random', 'chr22_KI270737v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270743v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1']
Keeping 81 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270736v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270538v1', 'chrUn_KI270589v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270744v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 81 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_plus_exact_body_0-mer.bw'
Merging 81 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:31. Running peak memory: 7.355GB.  
  PID: 25022;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.735GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam` (28234)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 7.355GB.  
  PID: 28234;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_minus_smooth_body_0-mer.bw -p 8 --variable-step --tail-edge --scale 21096561.0` (28252)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/aligned_hg38/H9_treated_PRO-seq_3_minus.bam'
Temporary files will be stored in: 'tmp_H9_treated_PRO-seq_3_minus_cuttrace_5fa0aciy'
Processing with 4 cores...
Discarding 109 chunk(s) of reads: ['chrM', 'chr1_KI270707v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270757v1']
Keeping 86 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr2_KI270716v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270589v1', 'chrUn_KI270584v1', 'chrUn_KI270593v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270756v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_minus_exact_body_0-mer.bw'
Merging 86 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/H9_treated_PRO-seq_3/signal_hg38/H9_treated_PRO-seq_3_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:07:35. Running peak memory: 7.355GB.  
  PID: 28252;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.754GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  3:16:06
*  Total elapsed time (all runs):  5:45:59
*         Peak memory (this run):  7.3549 GB
*        Pipeline completed time: 2020-06-15 10:33:24
