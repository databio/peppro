### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_70 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_70pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 24 -M 24000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-aj38-16c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/
*  Pipeline started at:   (06-11 17:07:30) elapsed: 9.0 _TIME_

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
*              `cores`:  `24`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_70pct.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `24000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `K562_PRO-seq_70`
*              `scale`:  `True`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `SINGLE`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `0`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_70pct.fastq.gz

> `File_mb`	25439.68	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:07:31) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_PRO_70pct.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz` (454878)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 454878;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1.fastq`  

> `pigz -f -p 24 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1.fastq` (454879)
<pre>
</pre>
Command completed. Elapsed time: 0:13:28. Running peak memory: 0.003GB.  
  PID: 454879;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	347627030	PEPPRO	_RES_

> `Fastq_reads`	347627030	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz']

### FASTQ processing:  (06-11 17:47:09) elapsed: 2378.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1_processed.fastq`  

> `(cutadapt -j 24 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt` (861)
<pre>
</pre>
Command completed. Elapsed time: 0:21:41. Running peak memory: 10.124GB.  
  PID: 861;	Command: cutadapt;	Return code: 0;	Memory used: 10.124GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1_processed.fastq` (6020,6021)
<pre>
</pre>
Command completed. Elapsed time: 0:10:35. Running peak memory: 10.124GB.  
  PID: 6020;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 6021;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	296156505.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	8705410.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.5042	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastqc/K562_PRO-seq_70_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (15005)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 10.124GB.  
  PID: 15005;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	338921620	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1_processed.fastq` (16673)
<pre>
Started analysis of K562_PRO-seq_70_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_70_R1_processed.fastq
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:13) elapsed: 4084.0 _TIME_
Child process 16673 (fastqc) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/recover.lock.fastqc__K562_PRO-seq_70_R1_processed_fastqc.html
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/recover.lock.trimmed_fastqc

### Pipeline failed at:  (06-11 18:55:14) elapsed: 0.0 _TIME_

Total time: 1:47:52
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_70 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_70pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 24 -M 24000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-32c1
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/
*  Pipeline started at:   (06-11 19:07:16) elapsed: 4.0 _TIME_

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
*              `cores`:  `24`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_70pct.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `24000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_PRO-seq_70`
*              `scale`:  `True`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `SINGLE`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `0`
*          `verbosity`:  `None`

----------------------------------------

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_70pct.fastq.gz

> `File_mb`	25439.68	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:07:17) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz']

### FASTQ processing:  (06-11 19:07:17) elapsed: 0.0 _TIME_


> `cutadapt --version`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1_processed.fastq`  
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/lock.fastqc__K562_PRO-seq_70_R1_processed_fastqc.html
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastqc/K562_PRO-seq_70_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (97835)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 97835;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	338921620	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/lock.trimmed_fastqc
Overwriting target...
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1_processed.fastq` (98154)
<pre>
Started analysis of K562_PRO-seq_70_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_70_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_70_R1_processed.fastq
Analysis complete for K562_PRO-seq_70_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:18:42. Running peak memory: 0.26GB.  
  PID: 98154;	Command: fastqc;	Return code: 0;	Memory used: 0.26GB

> `FastQC report r1`	fastqc/K562_PRO-seq_70_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/processed_R1.flag` (102892)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.26GB.  
  PID: 102892;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-11 19:28:38) elapsed: 1282.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt` (102899)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 0.317GB.  
  PID: 102899;	Command: Rscript;	Return code: 0;	Memory used: 0.317GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_70_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_70_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 19:28:43) elapsed: 4.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2313	PEPPRO	_RES_

### Prealignments (06-11 19:28:43) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:28:43) elapsed: 0.0 _TIME_


> `(bowtie2 -p 24 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq_70 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/prealignments/K562_PRO-seq_70_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
338921620 reads; of these:
  338921620 (100.00%) were unpaired; of these:
    307644743 (90.77%) aligned 0 times
    31276877 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	31276877.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-11 20:21:31) elapsed: 3168.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam`  

> `bowtie2 -p 24 --very-sensitive --rg-id K562_PRO-seq_70 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/prealignments/K562_PRO-seq_70_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/tmpt6y9en4j -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_temp.bam` (114135,114142,114143)
<pre>
307644743 reads; of these:
  307644743 (100.00%) were unpaired; of these:
    3774411 (1.23%) aligned 0 times
    229720636 (74.67%) aligned exactly 1 time
    74149696 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 99 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 2:28:26. Running peak memory: 4.018GB.  
  PID: 114142;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 114135;	Command: bowtie2;	Return code: 0;	Memory used: 4.018GB  
  PID: 114143;	Command: samtools;	Return code: 0;	Memory used: 0.904GB


> `samtools view -q 10 -b -@ 24 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam` (150154)
<pre>
</pre>
Command completed. Elapsed time: 0:07:25. Running peak memory: 4.018GB.  
  PID: 150154;	Command: samtools;	Return code: 0;	Memory used: 0.029GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	303870332	PEPPRO	_RES_

> `QC_filtered_reads`	32723809	PEPPRO	_RES_

> `Aligned_reads`	271146523	PEPPRO	_RES_

> `Alignment_rate`	80.0	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	21.5	PEPPRO	_RES_

### Compress all unmapped read files (06-11 23:42:36) elapsed: 12065.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/prealignments/K562_PRO-seq_70_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 24 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/prealignments/K562_PRO-seq_70_human_rDNA_unmap.fq` (155643)
<pre>
</pre>
Command completed. Elapsed time: 0:05:38. Running peak memory: 4.018GB.  
  PID: 155643;	Command: pigz;	Return code: 0;	Memory used: 0.021GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_temp.bam` (156194)
<pre>
</pre>
Command completed. Elapsed time: 0:04:37. Running peak memory: 4.018GB.  
  PID: 156194;	Command: samtools;	Return code: 0;	Memory used: 0.023GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	6405734	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam` (156679)
<pre>
</pre>
Command completed. Elapsed time: 0:04:06. Running peak memory: 4.018GB.  
  PID: 156679;	Command: samtools;	Return code: 0;	Memory used: 0.023GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/chr_sizes.bed` (157092,157093,157094,157095)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.018GB.  
  PID: 157094;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 157092;	Command: samtools;	Return code: 0;	Memory used: 0.009GB  
  PID: 157095;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 157093;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/chr_sizes.bed -b -@ 24 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_noMT.bam` (157097)
<pre>
</pre>
Command completed. Elapsed time: 0:04:15. Running peak memory: 4.018GB.  
  PID: 157097;	Command: samtools;	Return code: 0;	Memory used: 0.031GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam` (157605)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 4.018GB.  
  PID: 157605;	Command: mv;	Return code: 0;	Memory used: 0.001GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam` (157608)
<pre>
</pre>
Command completed. Elapsed time: 0:04:00. Running peak memory: 4.018GB.  
  PID: 157608;	Command: samtools;	Return code: 0;	Memory used: 0.024GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-12 00:17:06) elapsed: 2070.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam -c 24 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_bamQC.tsv` (159096)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/tmp_K562_PRO-seq_70_sort_43td5ke0'
Processing with 24 cores...
Discarding 83 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 112 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:04:26. Running peak memory: 6.149GB.  
  PID: 159096;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 6.149GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_bamQC.tsv`

> `NRF`	0.38	PEPPRO	_RES_

> `PBC1`	0.63	PEPPRO	_RES_

> `PBC2`	4.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_unmap.bam`  

> `samtools view -b -@ 24 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_unmap.bam` (159616)
<pre>
</pre>
Command completed. Elapsed time: 0:00:51. Running peak memory: 6.149GB.  
  PID: 159616;	Command: samtools;	Return code: 0;	Memory used: 0.014GB


> `samtools view -c -f 4 -@ 24 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_temp.bam`

> `Unmapped_reads`	3774411	PEPPRO	_RES_

### Split BAM by strand (06-12 00:23:12) elapsed: 366.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam` (159759)
<pre>
</pre>
Command completed. Elapsed time: 0:19:49. Running peak memory: 6.149GB.  
  PID: 159759;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam` (161682)
<pre>
</pre>
Command completed. Elapsed time: 0:17:54. Running peak memory: 6.149GB.  
  PID: 161682;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-12 01:00:55) elapsed: 2263.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (163482)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.149GB.  
  PID: 163482;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/plus_TSS.tsv -p ends -c 24 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_plus_TssEnrichment.txt` (163488)
<pre>
</pre>
Command completed. Elapsed time: 0:00:27. Running peak memory: 6.149GB.  
  PID: 163488;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.479GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/minus_TSS.tsv -p ends -c 24 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_minus_TssEnrichment.txt` (163579)
<pre>
</pre>
Command completed. Elapsed time: 0:00:22. Running peak memory: 6.149GB.  
  PID: 163579;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.755GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_minus_TssEnrichment.txt` (163648)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 6.149GB.  
  PID: 163648;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_70_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_70_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt` (163678,163679,163680,163681)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.149GB.  
  PID: 163678;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 163680;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 163679;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 163681;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_keep.txt` (163683)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.149GB.  
  PID: 163683;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-12 01:01:51) elapsed: 57.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_ensembl_tss.bed` (163685,163686)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 6.149GB.  
  PID: 163685;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 163686;	Command: bedtools;	Return code: 0;	Memory used: 0.049GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_ensembl_gene_body.bed` (163689,163690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 6.149GB.  
  PID: 163689;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 163690;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_TSS_density.bed` (163693,163694,163695,163696)
<pre>
</pre>
Command completed. Elapsed time: 0:05:49. Running peak memory: 6.149GB.  
  PID: 163696;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 163693;	Command: bedtools;	Return code: 0;	Memory used: 0.074GB  
  PID: 163695;	Command: sort;	Return code: 0;	Memory used: 0.012GB  
  PID: 163694;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_gene_body_density.bed` (164405,164406,164407)
<pre>
</pre>
Command completed. Elapsed time: 0:09:31. Running peak memory: 6.149GB.  
  PID: 164407;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 164405;	Command: bedtools;	Return code: 0;	Memory used: 0.722GB  
  PID: 164406;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/tmpgpevj5hw` (165608,165614,165615)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 6.149GB.  
  PID: 165608;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 165615;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 165614;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/tmpgpevj5hw | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-12 01:17:14) elapsed: 923.0 _TIME_

Total time: 6:10:03
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/tmpgpevj5hw | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_70 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_70pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 24 -M 24000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-14
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/
*  Pipeline started at:   (06-14 21:10:10) elapsed: 0.0 _TIME_

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
*              `cores`:  `24`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `True`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_70pct.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `24000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `True`
*        `sample_name`:  `K562_PRO-seq_70`
*              `scale`:  `True`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `SINGLE`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `0`
*          `verbosity`:  `None`

----------------------------------------

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_70pct.fastq.gz

> `File_mb`	25439.68	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:10:11) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/K562_PRO-seq_70_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/K562_PRO-seq_70.fastq.gz']

### FASTQ processing:  (06-14 21:10:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:10:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/cutadapt/K562_PRO-seq_70_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_70_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_70_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:10:11) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:10:11) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/prealignments/K562_PRO-seq_70_human_rDNA_unmap.fq

### Map to genome (06-14 21:10:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam`  

### Compress all unmapped read files (06-14 21:10:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/prealignments/K562_PRO-seq_70_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:10:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_unmap.bam`  

### Split BAM by strand (06-14 21:10:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam`  

### Calculate TSS enrichment (06-14 21:10:11) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_PRO-seq_70_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_70_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:10:11) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/tmp81sy5fh6` (110611,110612,110613)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 110611;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 110613;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 110612;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/tmp81sy5fh6 | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.156365) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/tmp81sy5fh6 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_pause_index.bed` (110619)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.006GB.  
  PID: 110619;	Command: awk;	Return code: 0;	Memory used: 0.004GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.88	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_pause_index.bed` (110624)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.204GB.  
  PID: 110624;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Pause index`	QC_hg38/K562_PRO-seq_70_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_70_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_pause_index.bed.gz`  

> `pigz -f -p 24 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_pause_index.bed` (110658)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.204GB.  
  PID: 110658;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:10:18) elapsed: 7.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam`
271146523 96714641

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam`
271146523 91984430

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_gene_sort.bed` (155332,155333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.204GB.  
  PID: 155333;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 155332;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_gene_coverage.bed` (155336)
<pre>
</pre>
Command completed. Elapsed time: 0:11:12. Running peak memory: 0.723GB.  
  PID: 155336;	Command: bedtools;	Return code: 0;	Memory used: 0.723GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/hg38_annotations.bed.gz` (5991)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.723GB.  
  PID: 5991;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 24 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/hg38_annotations.bed` (6036)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.723GB.  
  PID: 6036;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:30:05) elapsed: 1187.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/raw/hg38_annotations.bed` (6135)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.723GB.  
  PID: 6135;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Enhancer_sort.bed` (6407,6408,6409,6410)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.723GB.  
  PID: 6407;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 6408;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 6410;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 6409;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Enhancer_plus_coverage.bed` (6609)
<pre>
</pre>
Command completed. Elapsed time: 0:03:14. Running peak memory: 0.723GB.  
  PID: 6609;	Command: bedtools;	Return code: 0;	Memory used: 0.028GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Enhancer_minus_coverage.bed` (32302)
<pre>
</pre>
Command completed. Elapsed time: 0:03:03. Running peak memory: 0.723GB.  
  PID: 32302;	Command: bedtools;	Return code: 0;	Memory used: 0.05GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_sort.bed` (32655,32656,32657,32658)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.723GB.  
  PID: 32655;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 32656;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 32658;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB  
  PID: 32657;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_plus_coverage.bed` (32660)
<pre>
</pre>
Command completed. Elapsed time: 0:03:19. Running peak memory: 0.723GB.  
  PID: 32660;	Command: bedtools;	Return code: 0;	Memory used: 0.406GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_minus_coverage.bed` (32841)
<pre>
</pre>
Command completed. Elapsed time: 0:03:20. Running peak memory: 0.723GB.  
  PID: 32841;	Command: bedtools;	Return code: 0;	Memory used: 0.332GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_Flanking_Region"` (33244)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.723GB.  
  PID: 33244;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed` (33245,33246,33247,33248)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.723GB.  
  PID: 33245;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 33247;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 33246;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 33248;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_Flanking_Region_plus_coverage.bed` (33252)
<pre>
</pre>
Command completed. Elapsed time: 0:03:11. Running peak memory: 0.723GB.  
  PID: 33252;	Command: bedtools;	Return code: 0;	Memory used: 0.071GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_Flanking_Region_minus_coverage.bed` (33611)
<pre>
</pre>
Command completed. Elapsed time: 0:03:07. Running peak memory: 0.723GB.  
  PID: 33611;	Command: bedtools;	Return code: 0;	Memory used: 0.118GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/5_UTR"` (33773)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.723GB.  
  PID: 33773;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/5_UTR_sort.bed` (33775,33776,33777,33778)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.723GB.  
  PID: 33775;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 33776;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 33778;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 33777;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_5_UTR_plus_coverage.bed` (33781)
<pre>
</pre>
Command completed. Elapsed time: 0:03:05. Running peak memory: 0.723GB.  
  PID: 33781;	Command: bedtools;	Return code: 0;	Memory used: 0.056GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_5_UTR_minus_coverage.bed` (34174)
<pre>
</pre>
Command completed. Elapsed time: 0:03:01. Running peak memory: 0.723GB.  
  PID: 34174;	Command: bedtools;	Return code: 0;	Memory used: 0.046GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/3_UTR"` (34523)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.723GB.  
  PID: 34523;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/3_UTR_sort.bed` (34524,34525,34526,34527)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.723GB.  
  PID: 34524;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 34525;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 34527;	Command: bedtools;	Return code: 0;	Memory used: 0.034GB  
  PID: 34526;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_3_UTR_plus_coverage.bed` (34530)
<pre>
</pre>
Command completed. Elapsed time: 0:03:09. Running peak memory: 0.723GB.  
  PID: 34530;	Command: bedtools;	Return code: 0;	Memory used: 0.076GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_3_UTR_minus_coverage.bed` (34698)
<pre>
</pre>
Command completed. Elapsed time: 0:03:04. Running peak memory: 0.723GB.  
  PID: 34698;	Command: bedtools;	Return code: 0;	Memory used: 0.116GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Exon_sort.bed` (35150,35151,35152,35153)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.723GB.  
  PID: 35150;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 35151;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 35153;	Command: bedtools;	Return code: 0;	Memory used: 0.174GB  
  PID: 35152;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Exon_plus_coverage.bed` (35159)
<pre>
</pre>
Command completed. Elapsed time: 0:03:23. Running peak memory: 0.723GB.  
  PID: 35159;	Command: bedtools;	Return code: 0;	Memory used: 0.33GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Exon_minus_coverage.bed` (35523)
<pre>
</pre>
Command completed. Elapsed time: 0:03:21. Running peak memory: 0.723GB.  
  PID: 35523;	Command: bedtools;	Return code: 0;	Memory used: 0.133GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Intron_sort.bed` (35702,35703,35704,35705)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.723GB.  
  PID: 35702;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 35704;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 35703;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 35705;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Intron_plus_coverage.bed` (35710)
<pre>
</pre>
Command completed. Elapsed time: 0:04:06. Running peak memory: 0.723GB.  
  PID: 35710;	Command: bedtools;	Return code: 0;	Memory used: 0.321GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Intron_minus_coverage.bed` (58217)
<pre>
</pre>
Command completed. Elapsed time: 0:03:57. Running peak memory: 0.723GB.  
  PID: 58217;	Command: bedtools;	Return code: 0;	Memory used: 0.309GB


### Plot cFRiF/FRiF (06-14 22:16:36) elapsed: 2791.0 _TIME_


> `samtools view -@ 24 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_70 -z 3099922541 -n 134932628 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Intron_plus_coverage.bed` (58659)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 0.91GB.  
  PID: 58659;	Command: Rscript;	Return code: 0;	Memory used: 0.91GB

> `cFRiF`	QC_hg38/K562_PRO-seq_70_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_70_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_70 -z 3099922541 -n 134932628 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_Intron_plus_coverage.bed` (58711)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 0.91GB.  
  PID: 58711;	Command: Rscript;	Return code: 0;	Memory used: 0.432GB

> `FRiF`	QC_hg38/K562_PRO-seq_70_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_70_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 22:18:14) elapsed: 98.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_exons_sort.bed` (58748,58749)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.91GB.  
  PID: 58749;	Command: bedtools;	Return code: 0;	Memory used: 0.086GB  
  PID: 58748;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_introns_sort.bed` (58757,58758,58759)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.91GB.  
  PID: 58758;	Command: bedtools;	Return code: 0;	Memory used: 0.089GB  
  PID: 58757;	Command: grep;	Return code: 0;	Memory used: 0.004GB  
  PID: 58759;	Command: bedtools;	Return code: 0;	Memory used: 0.101GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exons_coverage.bed` (58766)
<pre>
</pre>
Command completed. Elapsed time: 0:06:40. Running peak memory: 0.91GB.  
  PID: 58766;	Command: bedtools;	Return code: 0;	Memory used: 0.207GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_introns_coverage.bed` (81559)
<pre>
</pre>
Command completed. Elapsed time: 0:08:34. Running peak memory: 0.91GB.  
  PID: 81559;	Command: bedtools;	Return code: 0;	Memory used: 0.345GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/271.146523)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exons_rpkm.bed` (112645,112712,112716)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.91GB.  
  PID: 112645;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 112716;	Command: sort;	Return code: 0;	Memory used: 0.002GB  
  PID: 112712;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/271.146523)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_introns_rpkm.bed` (113043,113045,113046)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.91GB.  
  PID: 113043;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 113046;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 113045;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exon_intron_ratios.bed` (113329,113330,113331)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.91GB.  
  PID: 113329;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 113331;	Command: sort;	Return code: 0;	Memory used: 0.006GB  
  PID: 113330;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exon_intron_ratios.bed --annotate` (113538)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.91GB.  
  PID: 113538;	Command: Rscript;	Return code: 0;	Memory used: 0.287GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_70_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_70_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exon_intron_ratios.bed.gz`  

> `pigz -f -p 24 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/QC_hg38/K562_PRO-seq_70_exon_intron_ratios.bed` (115401)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.91GB.  
  PID: 115401;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Produce bigWig files (06-14 22:33:48) elapsed: 934.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam` (115411)
<pre>
</pre>
Command completed. Elapsed time: 0:02:12. Running peak memory: 0.91GB.  
  PID: 115411;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_plus_smooth_body_0-mer.bw -p 16 --variable-step --tail-edge --scale 271146523.0` (126736)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_70_plus_cuttrace_bkj38iab'
Processing with 8 cores...
stdin is empty of data
Discarding 91 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1']
Keeping 104 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 104 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_plus_exact_body_0-mer.bw'
Merging 104 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:22. Running peak memory: 2.838GB.  
  PID: 126736;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.838GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam` (172665)
<pre>
</pre>
Command completed. Elapsed time: 0:02:21. Running peak memory: 2.838GB.  
  PID: 172665;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_minus_smooth_body_0-mer.bw -p 16 --variable-step --tail-edge --scale 271146523.0` (172828)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/aligned_hg38/K562_PRO-seq_70_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_70_minus_cuttrace_we5awi7m'
Processing with 8 cores...
stdin is empty of data
stdin is empty of data
Discarding 95 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 100 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 100 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_minus_exact_body_0-mer.bw'
Merging 100 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_70/signal_hg38/K562_PRO-seq_70_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:38. Running peak memory: 3.395GB.  
  PID: 172828;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.395GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:47:12
*  Total elapsed time (all runs):  10:19:25
*         Peak memory (this run):  3.3951 GB
*        Pipeline completed time: 2020-06-14 22:57:22
