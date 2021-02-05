### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_80 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_80pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 24 -M 24000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty`
*         Compute host:  udc-aj40-15c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/
*  Pipeline started at:   (06-11 17:07:30) elapsed: 8.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_80pct.fastq.gz']`
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
*        `sample_name`:  `K562_PRO-seq_80`
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

Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_80pct.fastq.gz

> `File_mb`	29053.1	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 17:07:31) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz`  

> `ln -sf /project/shefflab/data//guertin/fastq/K562_PRO_80pct.fastq.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz` (218552)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 218552;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1.fastq`  

> `pigz -f -p 24 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1.fastq` (218557)
<pre>
</pre>
Command completed. Elapsed time: 0:12:49. Running peak memory: 0.003GB.  
  PID: 218557;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	397272644	PEPPRO	_RES_

> `Fastq_reads`	397272644	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz']

### FASTQ processing:  (06-11 17:53:31) elapsed: 2760.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1_processed.fastq`  

> `(cutadapt -j 24 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1.fastq -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1_noadap.fastq --too-short-output /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1_short.fastq ) > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt` (228512)
<pre>
</pre>
Command completed. Elapsed time: 0:21:55. Running peak memory: 10.266GB.  
  PID: 228512;	Command: cutadapt;	Return code: 0;	Memory used: 10.266GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1_processed.fastq` (232399,232402)
<pre>
</pre>
Command completed. Elapsed time: 0:11:48. Running peak memory: 10.266GB.  
  PID: 232399;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 232402;	Command: seqtk;	Return code: 0;	Memory used: 0.002GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	338450269.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt`

> `wc -l /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1_short.fastq | awk '{print $1}'`

> `Uninformative_adapter_reads`	9949631.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	2.5045	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastqc/K562_PRO-seq_80_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (236482)
<pre>
### Calculate the number of trimmed reads
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 10.266GB.  
  PID: 236482;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	387323013	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1_processed.fastq` (236899)
<pre>
Started analysis of K562_PRO-seq_80_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_80_R1_processed.fastq
</pre>
Got SIGTERM. Failing gracefully... (06-11 18:55:14) elapsed: 3703.0 _TIME_
Child process 236899 (fastqc) terminated after 0 sec.
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/recover.lock.fastqc__K562_PRO-seq_80_R1_processed_fastqc.html
Setting dynamic recover file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/recover.lock.trimmed_fastqc

### Pipeline failed at:  (06-11 18:55:14) elapsed: 0.0 _TIME_

Total time: 1:47:52
Failure reason: SIGTERM
Pipeline aborted.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_80 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_80pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 24 -M 24000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba25-33c0
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/
*  Pipeline started at:   (06-11 19:07:16) elapsed: 6.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_80pct.fastq.gz']`
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
*        `sample_name`:  `K562_PRO-seq_80`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_80pct.fastq.gz

> `File_mb`	29053.1	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-11 19:07:17) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz']

### FASTQ processing:  (06-11 19:07:17) elapsed: 0.0 _TIME_


> `cutadapt --version`
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1_processed.fastq`  
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/lock.fastqc__K562_PRO-seq_80_R1_processed_fastqc.html
Overwriting target...
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastqc/K562_PRO-seq_80_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (4963)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0GB.  
  PID: 4963;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	387323013	PEPPRO	_RES_

> `Trim_loss_rate`	2.5	PEPPRO	_RES_
Found lock file: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/lock.trimmed_fastqc
Overwriting target...
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastqc /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1_processed.fastq` (5543)
<pre>
Started analysis of K562_PRO-seq_80_R1_processed.fastq
Approx 5% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 10% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 15% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 20% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 25% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 30% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 35% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 40% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 45% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 50% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 55% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 60% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 65% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 70% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 75% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 80% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 85% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 90% complete for K562_PRO-seq_80_R1_processed.fastq
Approx 95% complete for K562_PRO-seq_80_R1_processed.fastq
Analysis complete for K562_PRO-seq_80_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:21:17. Running peak memory: 0.277GB.  
  PID: 5543;	Command: fastqc;	Return code: 0;	Memory used: 0.277GB

> `FastQC report r1`	fastqc/K562_PRO-seq_80_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/processed_R1.flag` (8943)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.277GB.  
  PID: 8943;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (06-11 19:31:35) elapsed: 1459.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt` (8949)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.318GB.  
  PID: 8949;	Command: Rscript;	Return code: 0;	Memory used: 0.318GB

> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_80_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_80_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (06-11 19:31:40) elapsed: 5.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2313	PEPPRO	_RES_

### Prealignments (06-11 19:31:40) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-11 19:31:40) elapsed: 0.0 _TIME_


> `(bowtie2 -p 24 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /project/shefflab/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_PRO-seq_80 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/prealignments/K562_PRO-seq_80_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
387323013 reads; of these:
  387323013 (100.00%) were unpaired; of these:
    351574955 (90.77%) aligned 0 times
    35748058 (9.23%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
9.23% overall alignment rate

> `Aligned_reads_human_rDNA`	35748058.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	9.23	PEPPRO	_RES_

### Map to genome (06-11 20:29:46) elapsed: 3486.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam`  

> `bowtie2 -p 24 --very-sensitive --rg-id K562_PRO-seq_80 -x /project/shefflab/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/prealignments/K562_PRO-seq_80_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/tmpfrwarb4h -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_temp.bam` (18880,18887,18888)
<pre>
351574955 reads; of these:
  351574955 (100.00%) were unpaired; of these:
    4313427 (1.23%) aligned 0 times
    262525687 (74.67%) aligned exactly 1 time
    84735841 (24.10%) aligned >1 times
98.77% overall alignment rate
[bam_sort_core] merging from 113 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 2:48:39. Running peak memory: 4.028GB.  
  PID: 18880;	Command: bowtie2;	Return code: 0;	Memory used: 4.028GB  
  PID: 18887;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 18888;	Command: samtools;	Return code: 0;	Memory used: 0.909GB


> `samtools view -q 10 -b -@ 24 -U /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_fail_qc.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam` (55880)
<pre>
</pre>
Command completed. Elapsed time: 0:08:20. Running peak memory: 4.028GB.  
  PID: 55880;	Command: samtools;	Return code: 0;	Memory used: 0.028GB


> `samtools depth -b /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	347261528	PEPPRO	_RES_

> `QC_filtered_reads`	37392862	PEPPRO	_RES_

> `Aligned_reads`	309868666	PEPPRO	_RES_

> `Alignment_rate`	80.0	PEPPRO	_RES_

> `Total_efficiency`	78.0	PEPPRO	_RES_

> `Read_depth`	24.22	PEPPRO	_RES_

### Compress all unmapped read files (06-12 00:16:28) elapsed: 13602.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/prealignments/K562_PRO-seq_80_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 24 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/prealignments/K562_PRO-seq_80_human_rDNA_unmap.fq` (62213)
<pre>
</pre>
Command completed. Elapsed time: 0:06:13. Running peak memory: 4.028GB.  
  PID: 62213;	Command: pigz;	Return code: 0;	Memory used: 0.019GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_temp.bam` (62865)
<pre>
</pre>
Command completed. Elapsed time: 0:05:15. Running peak memory: 4.028GB.  
  PID: 62865;	Command: samtools;	Return code: 0;	Memory used: 0.026GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	7320857	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam` (63353)
<pre>
</pre>
Command completed. Elapsed time: 0:04:44. Running peak memory: 4.028GB.  
  PID: 63353;	Command: samtools;	Return code: 0;	Memory used: 0.022GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam | cut -f 1-2 | awk '{print $1, 0, $2}' | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/chr_sizes.bed` (63882,63883,63884,63885)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 4.028GB.  
  PID: 63884;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 63882;	Command: samtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 63885;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 63883;	Command: cut;	Return code: 0;	Memory used: 0.001GB


> `samtools view -L /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/chr_sizes.bed -b -@ 24 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_noMT.bam` (63888)
<pre>
</pre>
Command completed. Elapsed time: 0:04:35. Running peak memory: 4.028GB.  
  PID: 63888;	Command: samtools;	Return code: 0;	Memory used: 0.03GB


> `mv /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_noMT.bam /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam` (64390)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 4.028GB.  
  PID: 64390;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam` (64393)
<pre>
</pre>
Command completed. Elapsed time: 0:04:36. Running peak memory: 4.028GB.  
  PID: 64393;	Command: samtools;	Return code: 0;	Memory used: 0.022GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (06-12 00:55:39) elapsed: 2351.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio//peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam -c 24 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_bamQC.tsv` (66311)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/tmp_K562_PRO-seq_80_sort__6kt28qk'
Processing with 24 cores...
Discarding 81 chunk(s) of reads: ['chrM', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1']
Keeping 114 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
</pre>
Command completed. Elapsed time: 0:05:06. Running peak memory: 11.012GB.  
  PID: 66311;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamQC.py;	Return code: 0;	Memory used: 11.012GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_bamQC.tsv`

> `NRF`	0.34	PEPPRO	_RES_

> `PBC1`	0.6	PEPPRO	_RES_

> `PBC2`	4.14	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_unmap.bam`  

> `samtools view -b -@ 24 -f 4  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_temp.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_unmap.bam` (66917)
<pre>
</pre>
Command completed. Elapsed time: 0:00:59. Running peak memory: 11.012GB.  
  PID: 66917;	Command: samtools;	Return code: 0;	Memory used: 0.013GB


> `samtools view -c -f 4 -@ 24 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_temp.bam`

> `Unmapped_reads`	4313427	PEPPRO	_RES_

### Split BAM by strand (06-12 01:02:38) elapsed: 418.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam` (67086)
<pre>
</pre>
Command completed. Elapsed time: 0:22:26. Running peak memory: 11.012GB.  
  PID: 67086;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam` (71182)
<pre>
</pre>
Command completed. Elapsed time: 0:21:26. Running peak memory: 11.012GB.  
  PID: 71182;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (06-12 01:46:30) elapsed: 2632.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/minus_TSS.tsv' /project/shefflab/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (89156)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 11.012GB.  
  PID: 89156;	Command: sed;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/plus_TSS.tsv -p ends -c 24 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_plus_TssEnrichment.txt` (89162)
<pre>
</pre>
Command completed. Elapsed time: 0:00:32. Running peak memory: 11.012GB.  
  PID: 89162;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.291GB


> `TSS_coding_score`	13.9	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/minus_TSS.tsv -p ends -c 24 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_minus_TssEnrichment.txt` (89243)
<pre>
</pre>
Command completed. Elapsed time: 0:00:26. Running peak memory: 11.012GB.  
  PID: 89243;	Command: /scratch/jps3dp/tools/databio//peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 1.564GB


> `TSS_non-coding_score`	5.8	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_minus_TssEnrichment.txt` (89315)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_minus_TssEnrichment.txt
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
`geom_smooth()` using formula 'y ~ x'
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 11.012GB.  
  PID: 89315;	Command: Rscript;	Return code: 0;	Memory used: 0.316GB

> `TSS enrichment`	QC_hg38/K562_PRO-seq_80_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_80_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt` (89344,89345,89346,89347)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 11.012GB.  
  PID: 89344;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 89346;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 89345;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 89347;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_keep.txt` (89349)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 11.012GB.  
  PID: 89349;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (06-12 01:47:35) elapsed: 65.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_ensembl_tss.bed` (89351,89352)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 11.012GB.  
  PID: 89351;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 89352;	Command: bedtools;	Return code: 0;	Memory used: 0.051GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_ensembl_gene_body.bed` (89356,89357)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 11.012GB.  
  PID: 89356;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 89357;	Command: bedtools;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_TSS_density.bed` (89359,89360,89361,89362)
<pre>
</pre>
Command completed. Elapsed time: 0:06:36. Running peak memory: 11.012GB.  
  PID: 89361;	Command: sort;	Return code: 0;	Memory used: 0.014GB  
  PID: 89359;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB  
  PID: 89362;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 89360;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_gene_body_density.bed` (91773,91778,91779)
<pre>
</pre>
Command completed. Elapsed time: 0:10:54. Running peak memory: 11.012GB.  
  PID: 91773;	Command: bedtools;	Return code: 0;	Memory used: 0.823GB  
  PID: 91779;	Command: sort;	Return code: 0;	Memory used: 0.005GB  
  PID: 91778;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/tmpx39j3e7j` (95280,95287,95288)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 11.012GB.  
  PID: 95280;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 95288;	Command: env;	Return code: 0;	Memory used: 0.006GB  
  PID: 95287;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/tmpx39j3e7j | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}`
/bin/sh: -c: line 0: unexpected EOF while looking for matching `''
/bin/sh: -c: line 1: syntax error: unexpected end of file

### Pipeline failed at:  (06-12 02:05:07) elapsed: 1053.0 _TIME_

Total time: 6:57:58
Failure reason: Command 'awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/tmpx39j3e7j | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}' returned non-zero exit status 1.
### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio//peppro/pipelines/peppro.py --sample-name K562_PRO-seq_80 --genome hg38 --input /project/shefflab/data//guertin/fastq/K562_PRO_80pct.fastq.gz --single-or-paired SINGLE -O /project/shefflab/processed//peppro/paper/6.11.2020/results_pipeline -P 24 -M 24000 --protocol PRO --umi-len 0 --scale --prealignments human_rDNA --dirty -R`
*         Compute host:  udc-ba27-16
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/
*  Pipeline started at:   (06-14 21:10:20) elapsed: 8.0 _TIME_

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
*              `input`:  `['/project/shefflab/data//guertin/fastq/K562_PRO_80pct.fastq.gz']`
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
*        `sample_name`:  `K562_PRO-seq_80`
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

Removed existing flag: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/PEPPRO_failed.flag'
Local input file: /project/shefflab/data//guertin/fastq/K562_PRO_80pct.fastq.gz

> `File_mb`	29053.1	PEPPRO	_RES_

> `Read_type`	SINGLE	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (06-14 21:10:20) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz`  
Local input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz'
Found .fastq.gz file
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/K562_PRO-seq_80_R1.fastq`  
['/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/K562_PRO-seq_80.fastq.gz']

### FASTQ processing:  (06-14 21:10:20) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/fastq/processed_R1.flag`  

### Plot adapter insertion distribution (06-14 21:10:20) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/cutadapt/K562_PRO-seq_80_R1_adapter_insertion_distribution.pdf`  
> `Adapter insertion distribution`	cutadapt/K562_PRO-seq_80_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_PRO-seq_80_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_

### Prealignments (06-14 21:10:20) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (06-14 21:10:20) elapsed: 0.0 _TIME_

No files match cleanup pattern: /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/prealignments/K562_PRO-seq_80_human_rDNA_unmap.fq

### Map to genome (06-14 21:10:20) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam`  

### Compress all unmapped read files (06-14 21:10:20) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/prealignments/K562_PRO-seq_80_human_rDNA_unmap.fq.gz`  

### Calculate NRF, PBC1, and PBC2 (06-14 21:10:20) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam.bai`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_bamQC.tsv`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_unmap.bam`  

### Split BAM by strand (06-14 21:10:20) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam`  

### Calculate TSS enrichment (06-14 21:10:20) elapsed: 0.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_TSSenrichment.pdf`  
> `TSS enrichment`	QC_hg38/K562_PRO-seq_80_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_PRO-seq_80_TSSenrichment.png	PEPPRO	_OBJ_

### Calculate Pause Index (PI) (06-14 21:10:20) elapsed: 0.0 _TIME_

Missing stat 'Pause_index'
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_ensembl_tss.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_ensembl_gene_body.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_TSS_density.bed`  
Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_gene_body_density.bed`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_TSS_density.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_gene_body_density.bed | awk -v OFS='	' '{ if ($5 == "+"){print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($8-$2)^2))^2), ($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5} else {print $1, $2, $8, $4, sqrt((($6+$9)/sqrt(($3-$7)^2))^2),($6/sqrt(($3-$2)^2))/($9/sqrt(($8-$7)^2)), $5}}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/tmpulhj166k` (121324,121325,121326)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 121324;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 121326;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 121325;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/tmpulhj166k | sort -n | awk 'BEGIN{i=0} {s[i]=$1; i++;} END{print s[int(NR*0.5-0.5)]}'`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_pause_index.bed`  

> `awk -v OFS='	' '{ if ($5 > 0.17456) {print $1, $2, $3, $4, $6, $7}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/tmpulhj166k > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_pause_index.bed` (121333)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.004GB.  
  PID: 121333;	Command: awk;	Return code: 0;	Memory used: 0.004GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	10.91	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_pause_index.bed` (121338)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.285GB.  
  PID: 121338;	Command: Rscript;	Return code: 0;	Memory used: 0.285GB

> `Pause index`	QC_hg38/K562_PRO-seq_80_pause_index.pdf	Pause index	QC_hg38/K562_PRO-seq_80_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_pause_index.bed.gz`  

> `pigz -f -p 24 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_pause_index.bed` (121365)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.285GB.  
  PID: 121365;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (06-14 21:10:27) elapsed: 6.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam`
309868666 110527104

> `Plus_FRiP`	0.36	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam`
309868666 105122401

> `Minus_FRiP`	0.34	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_gene_sort.bed` (166353,166354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.285GB.  
  PID: 166354;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 166353;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_gene_coverage.bed` (166356)
<pre>
</pre>
Command completed. Elapsed time: 0:12:53. Running peak memory: 0.822GB.  
  PID: 166356;	Command: bedtools;	Return code: 0;	Memory used: 0.822GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/hg38_annotations.bed`  

> `ln -sf /project/shefflab/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/hg38_annotations.bed.gz` (37463)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.822GB.  
  PID: 37463;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 24 -d -c /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/hg38_annotations.bed` (37469)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.822GB.  
  PID: 37469;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (06-14 21:33:03) elapsed: 1357.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Enhancer`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/raw/hg38_annotations.bed` (37482)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 0.822GB.  
  PID: 37482;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Enhancer_sort.bed` (37485,37486,37487,37488)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.822GB.  
  PID: 37485;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 37486;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 37488;	Command: bedtools;	Return code: 0;	Memory used: 0.047GB  
  PID: 37487;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Enhancer_plus_coverage.bed` (37490)
<pre>
</pre>
Command completed. Elapsed time: 0:03:32. Running peak memory: 0.822GB.  
  PID: 37490;	Command: bedtools;	Return code: 0;	Memory used: 0.031GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Enhancer_minus_coverage.bed` (37864)
<pre>
</pre>
Command completed. Elapsed time: 0:03:30. Running peak memory: 0.822GB.  
  PID: 37864;	Command: bedtools;	Return code: 0;	Memory used: 0.056GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_sort.bed` (38272,38273,38274,38275)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.822GB.  
  PID: 38272;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 38273;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 38275;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 38274;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_plus_coverage.bed` (38278)
<pre>
</pre>
Command completed. Elapsed time: 0:03:48. Running peak memory: 0.822GB.  
  PID: 38278;	Command: bedtools;	Return code: 0;	Memory used: 0.475GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_minus_coverage.bed` (38481)
<pre>
</pre>
Command completed. Elapsed time: 0:03:47. Running peak memory: 0.822GB.  
  PID: 38481;	Command: bedtools;	Return code: 0;	Memory used: 0.38GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_Flanking_Region"` (38878)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.822GB.  
  PID: 38878;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_Flanking_Region_sort.bed` (38879,38880,38881,38882)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.822GB.  
  PID: 38879;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 38881;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 38880;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 38882;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_Flanking_Region_plus_coverage.bed` (38886)
<pre>
</pre>
Command completed. Elapsed time: 0:03:36. Running peak memory: 0.822GB.  
  PID: 38886;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_Flanking_Region_minus_coverage.bed` (39314)
<pre>
</pre>
Command completed. Elapsed time: 0:03:30. Running peak memory: 0.822GB.  
  PID: 39314;	Command: bedtools;	Return code: 0;	Memory used: 0.134GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/5_UTR"` (39497)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.822GB.  
  PID: 39497;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/5_UTR_sort.bed` (39498,39499,39500,39501)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.822GB.  
  PID: 39498;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 39499;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 39501;	Command: bedtools;	Return code: 0;	Memory used: 0.03GB  
  PID: 39500;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_5_UTR_plus_coverage.bed` (39505)
<pre>
</pre>
Command completed. Elapsed time: 0:03:30. Running peak memory: 0.822GB.  
  PID: 39505;	Command: bedtools;	Return code: 0;	Memory used: 0.063GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_5_UTR_minus_coverage.bed` (39879)
<pre>
</pre>
Command completed. Elapsed time: 0:03:23. Running peak memory: 0.822GB.  
  PID: 39879;	Command: bedtools;	Return code: 0;	Memory used: 0.04GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/3' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/3_UTR"` (40312)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.822GB.  
  PID: 40312;	Command: mv;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/3_UTR_sort.bed` (40313,40314,40315,40316)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.822GB.  
  PID: 40313;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 40314;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 40316;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 40315;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_3_UTR_plus_coverage.bed` (40318)
<pre>
</pre>
Command completed. Elapsed time: 0:03:38. Running peak memory: 0.822GB.  
  PID: 40318;	Command: bedtools;	Return code: 0;	Memory used: 0.085GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_3_UTR_minus_coverage.bed` (40746)
<pre>
</pre>
Command completed. Elapsed time: 0:03:28. Running peak memory: 0.822GB.  
  PID: 40746;	Command: bedtools;	Return code: 0;	Memory used: 0.131GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Exon_sort.bed` (40927,40928,40929,40930)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 0.822GB.  
  PID: 40927;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 40928;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 40930;	Command: bedtools;	Return code: 0;	Memory used: 0.174GB  
  PID: 40929;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Exon_plus_coverage.bed` (40934)
<pre>
</pre>
Command completed. Elapsed time: 0:04:03. Running peak memory: 0.822GB.  
  PID: 40934;	Command: bedtools;	Return code: 0;	Memory used: 0.376GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Exon_minus_coverage.bed` (63408)
<pre>
</pre>
Command completed. Elapsed time: 0:03:47. Running peak memory: 0.822GB.  
  PID: 63408;	Command: bedtools;	Return code: 0;	Memory used: 0.151GB

Target exists: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Intron_sort.bed` (63797,63798,63799,63800)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.822GB.  
  PID: 63797;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 63799;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 63798;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 63800;	Command: bedtools;	Return code: 0;	Memory used: 0.083GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Intron_plus_coverage.bed` (63803)
<pre>
</pre>
Command completed. Elapsed time: 0:04:33. Running peak memory: 0.822GB.  
  PID: 63803;	Command: bedtools;	Return code: 0;	Memory used: 0.363GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Intron_minus_coverage.bed` (64272)
<pre>
</pre>
Command completed. Elapsed time: 0:04:49. Running peak memory: 0.822GB.  
  PID: 64272;	Command: bedtools;	Return code: 0;	Memory used: 0.351GB


### Plot cFRiF/FRiF (06-14 22:26:07) elapsed: 3184.0 _TIME_


> `samtools view -@ 24 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_80 -z 3099922541 -n 154200778 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Intron_plus_coverage.bed` (86781)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:44. Running peak memory: 0.95GB.  
  PID: 86781;	Command: Rscript;	Return code: 0;	Memory used: 0.95GB

> `cFRiF`	QC_hg38/K562_PRO-seq_80_cFRiF.pdf	cFRiF	QC_hg38/K562_PRO-seq_80_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R frif -s K562_PRO-seq_80 -z 3099922541 -n 154200778 -y frif --reads -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Promoter_Flanking_Region_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_Intron_plus_coverage.bed` (86834)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:31. Running peak memory: 0.95GB.  
  PID: 86834;	Command: Rscript;	Return code: 0;	Memory used: 0.444GB

> `FRiF`	QC_hg38/K562_PRO-seq_80_FRiF.pdf	FRiF	QC_hg38/K562_PRO-seq_80_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (06-14 22:27:50) elapsed: 103.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_exons_sort.bed` (86871,86872)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.95GB.  
  PID: 86872;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB  
  PID: 86871;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `grep -wf /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_keep.txt /project/shefflab/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_introns_sort.bed` (86879,86880,86881)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 0.95GB.  
  PID: 86880;	Command: bedtools;	Return code: 0;	Memory used: 0.084GB  
  PID: 86879;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 86881;	Command: bedtools;	Return code: 0;	Memory used: 0.092GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exons_coverage.bed` (86891)
<pre>
</pre>
Command completed. Elapsed time: 0:08:06. Running peak memory: 0.95GB.  
  PID: 86891;	Command: bedtools;	Return code: 0;	Memory used: 0.236GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_sort.bam -g /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_introns_coverage.bed` (131852)
<pre>
</pre>
Command completed. Elapsed time: 0:10:31. Running peak memory: 0.95GB.  
  PID: 131852;	Command: bedtools;	Return code: 0;	Memory used: 0.392GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/309.868666)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exons_rpkm.bed` (182385,182389,182390)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.95GB.  
  PID: 182385;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 182390;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 182389;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/309.868666)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_introns_rpkm.bed` (182392,182393,182394)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 0.95GB.  
  PID: 182392;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 182394;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 182393;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_introns_rpkm.bed /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exon_intron_ratios.bed` (182397,182398,182399)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.95GB.  
  PID: 182397;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 182399;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 182398;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	1.37	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio//peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exon_intron_ratios.bed --annotate` (182405)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 0.95GB.  
  PID: 182405;	Command: Rscript;	Return code: 0;	Memory used: 0.319GB

> `mRNA contamination`	QC_hg38/K562_PRO-seq_80_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_PRO-seq_80_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exon_intron_ratios.bed.gz`  

> `pigz -f -p 24 -f /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/QC_hg38/K562_PRO-seq_80_exon_intron_ratios.bed` (182432)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.95GB.  
  PID: 182432;	Command: pigz;	Return code: 0;	Memory used: 0.005GB


### Produce bigWig files (06-14 22:46:47) elapsed: 1137.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam` (182441)
<pre>
</pre>
Command completed. Elapsed time: 0:02:23. Running peak memory: 0.95GB.  
  PID: 182441;	Command: samtools;	Return code: 0;	Memory used: 0.017GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_plus_smooth_body_0-mer.bw -p 16 --variable-step --tail-edge --scale 309868666.0` (5931)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_plus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_80_plus_cuttrace_exevh_2u'
Processing with 8 cores...
stdin is empty of data
Discarding 90 chunk(s) of reads: ['chrM', 'chr1_KI270710v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1']
Keeping 105 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270465v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270538v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270333v1', 'chrUn_KI270337v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 105 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_plus_exact_body_0-mer.bw'
Merging 105 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:09:52. Running peak memory: 3.114GB.  
  PID: 5931;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 3.114GB

Target to produce: `/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam` (54479)
<pre>
</pre>
Command completed. Elapsed time: 0:02:27. Running peak memory: 3.114GB.  
  PID: 54479;	Command: samtools;	Return code: 0;	Memory used: 0.016GB


> `/scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam -c /project/shefflab/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_minus_smooth_body_0-mer.bw -p 16 --variable-step --tail-edge --scale 309868666.0` (68770)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/aligned_hg38/K562_PRO-seq_80_minus.bam'
Temporary files will be stored in: 'tmp_K562_PRO-seq_80_minus_cuttrace_2s9ayk_y'
Processing with 8 cores...
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
stdin is empty of data
Discarding 94 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr9_KI270717v1_random', 'chr14_KI270723v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_KI270747v1', 'chrUn_KI270752v1', 'chrUn_KI270755v1']
Keeping 101 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr5_GL000208v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chrUn_KI270310v1', 'chrUn_KI270411v1', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270515v1', 'chrUn_KI270538v1', 'chrUn_KI270583v1', 'chrUn_KI270581v1', 'chrUn_KI270589v1', 'chrUn_KI270588v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_KI270336v1', 'chrUn_KI270448v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_GL000213v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270753v1', 'chrUn_KI270754v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1', 'chrEBV']
Reduce step (merge files)...
Merging 101 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_minus_exact_body_0-mer.bw'
Merging 101 files into output file: '/project/shefflab/processed/peppro/paper/6.11.2020/results_pipeline/K562_PRO-seq_80/signal_hg38/K562_PRO-seq_80_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:10:17. Running peak memory: 3.114GB.  
  PID: 68770;	Command: /scratch/jps3dp/tools/databio//peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 2.792GB


### Pipeline completed. Epilogue
*        Elapsed time (this run):  2:01:34
*  Total elapsed time (all runs):  11:38:48
*         Peak memory (this run):  3.1136 GB
*        Pipeline completed time: 2020-06-14 23:11:45
