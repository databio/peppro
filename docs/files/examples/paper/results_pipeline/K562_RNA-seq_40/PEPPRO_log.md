### Pipeline run code and environment:

*              Command:  `/scratch/jps3dp/tools/databio/peppro/pipelines/peppro.py --sample-name K562_RNA-seq_40 --genome hg38 --input /project/shefflab/data/guertin/fastq/K562_40pct_RNArc_r2.fastq.gz --single-or-paired single --protocol PRO --prealignments human_rDNA -O /project/shefflab/processed/peppro/paper/dev4/results_pipeline -P 4 -M 16000`
*         Compute host:  udc-ba26-10
*          Working dir:  /sfs/lustre/bahamut/scratch/jps3dp/tools/databio/ppqc
*            Outfolder:  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/
*  Pipeline started at:   (02-27 09:23:32) elapsed: 0.0 _TIME_

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
*              `cores`:  `4`
*           `coverage`:  `False`
*              `dedup`:  `seqkit`
*              `dirty`:  `False`
*  `ensembl_gene_body`:  `None`
*        `ensembl_tss`:  `None`
*          `exon_name`:  `None`
*       `force_follow`:  `False`
*    `genome_assembly`:  `hg38`
*              `input`:  `['/project/shefflab/data/guertin/fastq/K562_40pct_RNArc_r2.fastq.gz']`
*             `input2`:  `None`
*        `intron_name`:  `None`
*               `keep`:  `False`
*             `logdev`:  `False`
*            `max_len`:  `-1`
*                `mem`:  `16000`
*          `new_start`:  `False`
*            `no_fifo`:  `False`
*      `output_parent`:  `/project/shefflab/processed/peppro/paper/dev4/results_pipeline`
*         `paired_end`:  `False`
*           `pre_name`:  `None`
*      `prealignments`:  `['human_rDNA']`
*         `prioritize`:  `False`
*           `protocol`:  `PRO`
*            `recover`:  `False`
*        `sample_name`:  `K562_RNA-seq_40`
*              `scale`:  `False`
*        `search_file`:  `None`
*             `silent`:  `False`
*   `single_or_paired`:  `single`
*                `sob`:  `False`
*           `testmode`:  `False`
*            `trimmer`:  `seqtk`
*            `umi_len`:  `0`
*          `verbosity`:  `None`

----------------------------------------

Local input file: /project/shefflab/data/guertin/fastq/K562_40pct_RNArc_r2.fastq.gz

> `File_mb`	800.5	PEPPRO	_RES_

> `Read_type`	single	PEPPRO	_RES_

> `Genome`	hg38	PEPPRO	_RES_
Detected PRO input

### Merge/link and fastq conversion:  (02-27 09:23:33) elapsed: 1.0 _TIME_

Number of input file sets: 1
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz`  

> `ln -sf /project/shefflab/data/guertin/fastq/K562_40pct_RNArc_r2.fastq.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz` (142017)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.0GB.  
  PID: 142017;	Command: ln;	Return code: 0;	Memory used: 0.0GB

Local input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz'
Found .fastq.gz file
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1.fastq`  

> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1.fastq` (142018)
<pre>
</pre>
Command completed. Elapsed time: 0:00:20. Running peak memory: 0.003GB.  
  PID: 142018;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


> `Raw_reads`	10000000	PEPPRO	_RES_

> `Fastq_reads`	10000000	PEPPRO	_RES_
['/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/K562_RNA-seq_40.fastq.gz']

### FASTQ processing:  (02-27 09:24:04) elapsed: 31.0 _TIME_


> `cutadapt --version`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_processed.fastq`  

> `(cutadapt -j 4 -m 2 -O 1 -a TGGAATTCTCGGGTGCCAAGG /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1.fastq -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_noadap.fastq) > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt` (142064)
<pre>
</pre>
Command completed. Elapsed time: 0:00:34. Running peak memory: 0.163GB.  
  PID: 142064;	Command: cutadapt;	Return code: 0;	Memory used: 0.163GB


> `seqtk trimfq -b 0 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_noadap.fastq | seqtk seq -L 2 -r - > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_processed.fastq` (142108,142109)
<pre>
</pre>
Command completed. Elapsed time: 0:00:14. Running peak memory: 0.163GB.  
  PID: 142108;	Command: seqtk;	Return code: 0;	Memory used: 0.001GB  
  PID: 142109;	Command: seqtk;	Return code: 0;	Memory used: 0.003GB


> `grep 'Reads with adapters:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Reads_with_adapter`	6515223.0	PEPPRO	_RES_

> `grep 'Total basepairs processed:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `awk '{sum+=$1*$2} END {printf "%.0f", sum}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt`

> `grep 'Reads that were too short:' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{print $(NF-1)}'`

> `Uninformative_adapter_reads`	149598.0	PEPPRO	_RES_

> `Pct_uninformative_adapter_reads`	1.496	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastqc/K562_RNA-seq_40_R1_processed_fastqc.html`  

> `echo '### Calculate the number of trimmed reads'` (142143)
### Calculate the number of trimmed reads
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.163GB.  
  PID: 142143;	Command: echo;	Return code: 0;	Memory used: 0.0GB

Evaluating read trimming

> `Trimmed_reads`	9850402	PEPPRO	_RES_

> `Trim_loss_rate`	1.5	PEPPRO	_RES_
Targetless command, running...  

> `fastqc --noextract --outdir /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastqc /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_processed.fastq` (142149)
<pre>
Picked up JAVA_TOOL_OPTIONS: -Xss1280k
Picked up _JAVA_OPTIONS: -Djava.io.tmpdir=/scratch/jps3dp/tmp
Started analysis of K562_RNA-seq_40_R1_processed.fastq
Approx 5% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 10% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 15% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 20% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 25% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 30% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 35% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 40% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 45% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 50% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 55% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 60% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 65% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 70% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 75% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 80% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 85% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 90% complete for K562_RNA-seq_40_R1_processed.fastq
Approx 95% complete for K562_RNA-seq_40_R1_processed.fastq
Analysis complete for K562_RNA-seq_40_R1_processed.fastq
</pre>
Command completed. Elapsed time: 0:00:50. Running peak memory: 0.191GB.  
  PID: 142149;	Command: fastqc;	Return code: 0;	Memory used: 0.191GB

> `FastQC report r1`	fastqc/K562_RNA-seq_40_R1_processed_fastqc.html	FastQC report r1	None	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/processed_R1.flag`  

> `touch /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/processed_R1.flag` (142403)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 0.191GB.  
  PID: 142403;	Command: touch;	Return code: 0;	Memory used: 0.002GB


### Plot adapter insertion distribution (02-27 09:25:46) elapsed: 102.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_adapter_insertion_distribution.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R cutadapt -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt` (142405)
<pre>
Adapter insertion distribution plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 0.204GB.  
  PID: 142405;	Command: Rscript;	Return code: 0;	Memory used: 0.204GB

> `Adapter insertion distribution`	cutadapt/K562_RNA-seq_40_R1_adapter_insertion_distribution.pdf	Adapter insertion distribution	cutadapt/K562_RNA-seq_40_R1_adapter_insertion_distribution.png	PEPPRO	_OBJ_
Missing stat 'Peak_adapter_insertion_size'

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk 'BEGIN{max=   0; max_len=0; len=0}{if ($2>0+max) {max=$2; len=$1}; max_len=$1} END{print max_len-len}'`

> `Peak_adapter_insertion_size`	34	PEPPRO	_RES_
Missing stat 'Degradation_ratio'

###  Calculating degradation ratio (02-27 09:25:52) elapsed: 6.0 _TIME_


> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 10) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 20) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 30) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk 'NR>2 {print prev} {prev=$0}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk '{ if ($1 == 40) {status = 1}} END {if (status) {print status} else {print 0}}'`

> `awk '/count/,0' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/cutadapt/K562_RNA-seq_40_R1_cutadapt.txt | awk 'NR>2 {print prev} {prev=$0}' | awk '{if ($3/$2 < 0.01) print $1, $2}' | awk '{a[NR]=$1; b[NR]=$2; max_len=$1}{if ($1 > max_len) {max_len=$1}} END{ for (i in a) print 1+max_len-a[i], b[i]}' | sort -nk1 | awk '($1 <= 20 && $1 >= 10){degradedSum += $2}; ($1 >= 30 && $1 <= 40){intactSum += $2} END {if (intactSum < 1) {intactSum = 1} print degradedSum/intactSum}'`

> `Degradation_ratio`	0.2301	PEPPRO	_RES_

### Prealignments (02-27 09:25:52) elapsed: 0.0 _TIME_

Prealignment assemblies: ['human_rDNA']

### Map to human_rDNA (02-27 09:25:52) elapsed: 0.0 _TIME_


> `(bowtie2 -p 4 -k 1 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -x /scratch/jps3dp/DATA/genomes/human_rDNA/bowtie2_index/default/human_rDNA --rg-id K562_RNA-seq_40 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/fastq/K562_RNA-seq_40_R1_processed.fastq --un /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/prealignments/K562_RNA-seq_40_human_rDNA_unmap.fq 2>&1 > /dev/null)`
Missing stat 'Aligned_reads_human_rDNA'
9850402 reads; of these:
  9850402 (100.00%) were unpaired; of these:
    9212394 (93.52%) aligned 0 times
    638008 (6.48%) aligned exactly 1 time
    0 (0.00%) aligned >1 times
6.48% overall alignment rate

> `Aligned_reads_human_rDNA`	638008.0	PEPPRO	_RES_

> `Alignment_rate_human_rDNA`	6.48	PEPPRO	_RES_

### Map to genome (02-27 09:27:29) elapsed: 96.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam`  

> `bowtie2 -p 4 --very-sensitive -X 2000 --rg-id K562_RNA-seq_40 -x /scratch/jps3dp/DATA/genomes/hg38/bowtie2_index/default/hg38 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/prealignments/K562_RNA-seq_40_human_rDNA_unmap.fq | samtools view -bS - -@ 1  | samtools sort - -@ 1 -T /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/tmpjaiq_t7o -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam` (142540,142541,142542)
<pre>
9212394 reads; of these:
  9212394 (100.00%) were unpaired; of these:
    663304 (7.20%) aligned 0 times
    5915800 (64.22%) aligned exactly 1 time
    2633290 (28.58%) aligned >1 times
92.80% overall alignment rate
[bam_sort_core] merging from 2 files and 1 in-memory blocks...
</pre>
Command completed. Elapsed time: 0:16:03. Running peak memory: 3.534GB.  
  PID: 142540;	Command: bowtie2;	Return code: 0;	Memory used: 3.534GB  
  PID: 142541;	Command: samtools;	Return code: 0;	Memory used: 0.004GB  
  PID: 142542;	Command: samtools;	Return code: 0;	Memory used: 0.865GB


> `samtools view -q 10 -b -@ 4 -U /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_fail_qc.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam` (144052)
<pre>
</pre>
Command completed. Elapsed time: 0:00:37. Running peak memory: 3.534GB.  
  PID: 144052;	Command: samtools;	Return code: 0;	Memory used: 0.012GB


> `samtools depth -b /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam | awk '{counter++;sum+=$3}END{print sum/counter}'`

> `Mapped_reads`	8549090	PEPPRO	_RES_

> `QC_filtered_reads`	1305366	PEPPRO	_RES_

> `Aligned_reads`	7243724	PEPPRO	_RES_

> `Alignment_rate`	73.54	PEPPRO	_RES_

> `Total_efficiency`	72.44	PEPPRO	_RES_

> `Read_depth`	2.26	PEPPRO	_RES_

### Compress all unmapped read files (02-27 09:52:04) elapsed: 1475.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/prealignments/K562_RNA-seq_40_human_rDNA_unmap.fq.gz`  

> `pigz -f -p 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/prealignments/K562_RNA-seq_40_human_rDNA_unmap.fq` (144920)
<pre>
</pre>
Command completed. Elapsed time: 0:00:52. Running peak memory: 3.534GB.  
  PID: 144920;	Command: pigz;	Return code: 0;	Memory used: 0.004GB

Missing stat 'Mitochondrial_reads'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam.bai`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam` (144983)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.534GB.  
  PID: 144983;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam | grep -we 'chrM' -we 'chrMT' -we 'M' -we 'MT' -we 'rCRSd' -we 'rCRSd_3k'| cut -f 3`

> `Mitochondrial_reads`	351192	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_noMT.bam`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam` (144998)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 144998;	Command: samtools;	Return code: 0;	Memory used: 0.009GB


> `samtools idxstats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam | cut -f 1 | grep -vwe 'chrM' -vwe 'chrMT' -vwe 'M' -vwe 'MT' -vwe 'rCRSd' -vwe 'rCRSd_3k'| xargs samtools view -b -@ 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_noMT.bam` (145009,145010,145011,145012)
<pre>
</pre>
Command completed. Elapsed time: 0:00:21. Running peak memory: 3.534GB.  
  PID: 145010;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 145009;	Command: samtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 145011;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 145012;	Command: xargs;	Return code: 0;	Memory used: 0.027GB


> `mv /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_noMT.bam /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam` (145038)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 145038;	Command: mv;	Return code: 0;	Memory used: 0.002GB


> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam` (145039)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.534GB.  
  PID: 145039;	Command: samtools;	Return code: 0;	Memory used: 0.008GB

Missing stat 'Maximum_read_length'

> `samtools stats /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam | grep '^SN' | cut -f 2- | grep 'maximum length:' | cut -f 2-`

> `Maximum_read_length`	100	PEPPRO	_RES_
Missing stat 'Genome_size'

> `awk '{sum+=$2} END {printf "%.0f", sum}' /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes`

> `Genome_size`	3099922541	PEPPRO	_RES_

### Calculate NRF, PBC1, and PBC2 (02-27 09:54:12) elapsed: 128.0 _TIME_

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam.bai`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv`  

> `/scratch/jps3dp/tools/databio/peppro/tools/bamQC.py --silent -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -c 4 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv` (145076)
<pre>
Configured logger 'root' using pararead v0.6
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam'
Temporary files will be stored in: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/tmp_K562_RNA-seq_40_sort_4_ws3wu_'
Processing with 4 cores...
Discarding 117 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr22_KI270731v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrEBV']
Keeping 78 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270580v1', 'chrUn_KI270589v1', 'chrUn_KI270582v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
</pre>
Command completed. Elapsed time: 0:00:19. Running peak memory: 3.534GB.  
  PID: 145076;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamQC.py;	Return code: 0;	Memory used: 0.401GB


> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "NRF") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC1") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv`

> `awk '{ for (i=1; i<=NF; ++i) { if ($i ~ "PBC2") c=i } getline; print $c }' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_bamQC.tsv`

> `NRF`	0.91	PEPPRO	_RES_

> `PBC1`	0.97	PEPPRO	_RES_

> `PBC2`	41.7	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_unmap.bam`  

> `samtools view -b -@ 4 -f 4  /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_unmap.bam` (145109)
<pre>
</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 145109;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `samtools view -c -f 4 -@ 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_temp.bam`

> `Unmapped_reads`	663304	PEPPRO	_RES_

### Split BAM by strand (02-27 09:54:39) elapsed: 27.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam`  

> `samtools view -bh -F 20 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam` (145129)
<pre>
</pre>
Command completed. Elapsed time: 0:00:41. Running peak memory: 3.534GB.  
  PID: 145129;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


> `samtools view -bh -f 16 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam` (145354)
<pre>
</pre>
Command completed. Elapsed time: 0:00:40. Running peak memory: 3.534GB.  
  PID: 145354;	Command: samtools;	Return code: 0;	Memory used: 0.006GB


### Calculate TSS enrichment (02-27 09:56:00) elapsed: 81.0 _TIME_

Missing stat 'TSS_non-coding_score'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/plus_TSS.tsv`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/minus_TSS.tsv`  

> `sed -n -e '/[[:space:]]+/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/plus_TSS.tsv' -e '/[[:space:]]-/w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/minus_TSS.tsv' /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_TSS.bed` (145404)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 145404;	Command: sed;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_plus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/plus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_plus_TssEnrichment.txt` (145406)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 145406;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.196GB


> `TSS_coding_score`	16.7	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_minus_TssEnrichment.txt`  

> `/scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/minus_TSS.tsv -p ends -c 4 -z -v -s 6 -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_minus_TssEnrichment.txt` (145423)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 145423;	Command: /scratch/jps3dp/tools/databio/peppro/tools/pyTssEnrichment.py;	Return code: 0;	Memory used: 0.2GB


> `TSS_non-coding_score`	5.0	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSSenrichment.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R tss -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_plus_TssEnrichment.txt /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_minus_TssEnrichment.txt` (145443)
<pre>

Generating TSS plot with /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_plus_TssEnrichment.txt and /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_minus_TssEnrichment.txt
TSS enrichment plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 145443;	Command: Rscript;	Return code: 0;	Memory used: 0.221GB

> `TSS enrichment`	QC_hg38/K562_RNA-seq_40_TSSenrichment.pdf	TSS enrichment	QC_hg38/K562_RNA-seq_40_TSSenrichment.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt`  

> `samtools view -H /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam | grep 'SN:' | awk -F':' '{print $2,$3}' | awk -F' ' -v OFS='	' '{print $1,$3}' > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt` (145461,145462,145463,145464)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 145461;	Command: samtools;	Return code: 0;	Memory used: 0.0GB  
  PID: 145463;	Command: awk;	Return code: 0;	Memory used: 0.0GB  
  PID: 145462;	Command: grep;	Return code: 0;	Memory used: 0.0GB  
  PID: 145464;	Command: awk;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt` (145466)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 145466;	Command: cut;	Return code: 0;	Memory used: 0.002GB


### Calculate Pause Index (PI) (02-27 09:56:22) elapsed: 22.0 _TIME_

Missing stat 'Pause_index'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_tss.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_gene_body.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_TSS.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_tss.bed` (145468,145469)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.534GB.  
  PID: 145468;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 145469;	Command: bedtools;	Return code: 0;	Memory used: 0.093GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/ensembl_gtf/default/hg38_ensembl_gene_body.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_gene_body.bed` (145473,145474)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.534GB.  
  PID: 145473;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 145474;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSS_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_tss.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4,4 -k7,7nr | sort -k4,4 -u > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSS_density.bed` (145478,145479,145480,145481)
<pre>
</pre>
Command completed. Elapsed time: 0:00:16. Running peak memory: 3.534GB.  
  PID: 145478;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB  
  PID: 145480;	Command: sort;	Return code: 0;	Memory used: 0.007GB  
  PID: 145479;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 145481;	Command: sort;	Return code: 0;	Memory used: 0.011GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_gene_body_density.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_ensembl_gene_body.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | awk '$7>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_gene_body_density.bed` (145498,145499,145500)
<pre>
</pre>
Command completed. Elapsed time: 0:00:30. Running peak memory: 3.534GB.  
  PID: 145499;	Command: awk;	Return code: 0;	Memory used: 0.001GB  
  PID: 145498;	Command: bedtools;	Return code: 0;	Memory used: 0.036GB  
  PID: 145500;	Command: sort;	Return code: 0;	Memory used: 0.005GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed`  

> `join --nocheck-order -j4 -o 1.1 1.2 1.3 1.4 1.6 1.7 2.2 2.3 2.7 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_TSS_density.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_gene_body_density.bed | awk -v OFS='	' '{print $1, $2, $3, $4, ($6/($3-$2))/($9/($8-$7)), $5}' | env LC_COLLATE=C sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed` (145531,145532,145533)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 145531;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 145533;	Command: env;	Return code: 0;	Memory used: 0.004GB  
  PID: 145532;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `sort -k5,5n /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed | awk ' { a[i++]=$5; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `Pause_index`	12.91	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R pi --annotate -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed` (145538)
<pre>
Pause index plot completed!

</pre>
Command completed. Elapsed time: 0:00:05. Running peak memory: 3.534GB.  
  PID: 145538;	Command: Rscript;	Return code: 0;	Memory used: 0.215GB

> `Pause index`	QC_hg38/K562_RNA-seq_40_pause_index.pdf	Pause index	QC_hg38/K562_RNA-seq_40_pause_index.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_pause_index.bed` (145558)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 145558;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Calculate Fraction of Reads in pre-mature mRNA (02-27 09:57:21) elapsed: 58.0 _TIME_

Missing stat 'Plus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam`
7243724 2698895

> `Plus_FRiP`	0.37	PEPPRO	_RES_
Missing stat 'Minus_FRiP'

> `samtools view -@ 4 -c -L /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam`
7243724 2633125

> `Minus_FRiP`	0.36	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_gene_coverage.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_pre-mRNA.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_gene_sort.bed` (145588,145589)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 145589;	Command: bedtools;	Return code: 0;	Memory used: 0.013GB  
  PID: 145588;	Command: grep;	Return code: 0;	Memory used: 0.004GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_gene_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_gene_coverage.bed` (145592)
<pre>
</pre>
Command completed. Elapsed time: 0:00:25. Running peak memory: 3.534GB.  
  PID: 145592;	Command: bedtools;	Return code: 0;	Memory used: 0.039GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed`  

> `ln -sf /scratch/jps3dp/DATA/genomes/hg38/feat_annotation/default/hg38_annotations.bed.gz /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed.gz` (145615)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 145615;	Command: ln;	Return code: 0;	Memory used: 0.002GB


> `pigz -f -p 4 -d -c /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed.gz > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed` (145616)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 145616;	Command: pigz;	Return code: 0;	Memory used: 0.003GB


### Calculate cumulative and terminal fraction of reads in features (cFRiF/FRiF) (02-27 09:58:02) elapsed: 41.0 _TIME_


> `cut -f 4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed | sort -u`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/3' UTR`  

> `awk -F'	' '{print>"/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/"$4}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/raw/hg38_annotations.bed` (145625)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.534GB.  
  PID: 145625;	Command: awk;	Return code: 0;	Memory used: 0.002GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/3' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR"` (145629)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 145629;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR_sort.bed` (145630,145631,145632,145633)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 145630;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 145631;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 145633;	Command: bedtools;	Return code: 0;	Memory used: 0.035GB  
  PID: 145632;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_plus_coverage.bed` (145635)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 145635;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/3_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_minus_coverage.bed` (145644)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 145644;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/5' UTR`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/5' UTR" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR"` (145652)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 145652;	Command: mv;	Return code: 0;	Memory used: 0.0GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR_sort.bed` (145653,145654,145655,145656)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 145653;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 145654;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 145656;	Command: bedtools;	Return code: 0;	Memory used: 0.007GB  
  PID: 145655;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_plus_coverage.bed` (145659)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 145659;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/5_UTR_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_minus_coverage.bed` (145669)
<pre>
</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 145669;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer_sort.bed` (145676,145677,145679,145680)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 145676;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 145677;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 145680;	Command: bedtools;	Return code: 0;	Memory used: 0.045GB  
  PID: 145679;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_plus_coverage.bed` (145682)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 145682;	Command: bedtools;	Return code: 0;	Memory used: 0.009GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Enhancer_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_minus_coverage.bed` (145690)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 145690;	Command: bedtools;	Return code: 0;	Memory used: 0.01GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon_sort.bed` (145698,145699,145700,145701)
<pre>
</pre>
Command completed. Elapsed time: 0:00:03. Running peak memory: 3.534GB.  
  PID: 145698;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 145699;	Command: grep;	Return code: 0;	Memory used: 0.003GB  
  PID: 145701;	Command: bedtools;	Return code: 0;	Memory used: 0.161GB  
  PID: 145700;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_plus_coverage.bed` (145706)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.534GB.  
  PID: 145706;	Command: bedtools;	Return code: 0;	Memory used: 0.021GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Exon_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_minus_coverage.bed` (145716)
<pre>
</pre>
Command completed. Elapsed time: 0:00:10. Running peak memory: 3.534GB.  
  PID: 145716;	Command: bedtools;	Return code: 0;	Memory used: 0.014GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron_sort.bed` (145729,145731,145732,145733)
<pre>
</pre>
Command completed. Elapsed time: 0:00:02. Running peak memory: 3.534GB.  
  PID: 145729;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 145732;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 145731;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 145733;	Command: bedtools;	Return code: 0;	Memory used: 0.079GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_plus_coverage.bed` (145736)
<pre>
</pre>
Command completed. Elapsed time: 0:00:11. Running peak memory: 3.534GB.  
  PID: 145736;	Command: bedtools;	Return code: 0;	Memory used: 0.024GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Intron_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_minus_coverage.bed` (145747)
<pre>
</pre>
Command completed. Elapsed time: 0:00:09. Running peak memory: 3.534GB.  
  PID: 145747;	Command: bedtools;	Return code: 0;	Memory used: 0.019GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_sort.bed` (145756,145757,145758,145759)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 145756;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 145757;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 145759;	Command: bedtools;	Return code: 0;	Memory used: 0.006GB  
  PID: 145758;	Command: cut;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_plus_coverage.bed` (145761)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 145761;	Command: bedtools;	Return code: 0;	Memory used: 0.02GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_minus_coverage.bed` (145770)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 145770;	Command: bedtools;	Return code: 0;	Memory used: 0.015GB

Target exists: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter Flanking Region`  
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region`  

> `mv "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter Flanking Region" "/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region"` (145778)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 145778;	Command: mv;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed`  

> `cut -f 1 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | grep -wf - /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region | cut -f 1-3 | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed` (145779,145780,145781,145782)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 145779;	Command: cut;	Return code: 0;	Memory used: 0.0GB  
  PID: 145781;	Command: cut;	Return code: 0;	Memory used: 0.001GB  
  PID: 145780;	Command: grep;	Return code: 0;	Memory used: 0.002GB  
  PID: 145782;	Command: bedtools;	Return code: 0;	Memory used: 0.048GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_plus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_plus_coverage.bed` (145785)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 145785;	Command: bedtools;	Return code: 0;	Memory used: 0.008GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_minus_coverage.bed`  

> `bedtools coverage -sorted  -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/Promoter_Flanking_Region_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_minus_coverage.bed` (145903)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 145903;	Command: bedtools;	Return code: 0;	Memory used: 0.012GB


### Plot cFRiF/FRiF (02-27 10:00:09) elapsed: 127.0 _TIME_


> `samtools view -@ 4 -q 10 -c -F4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam`
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_cFRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_40 -z 3099922541 -n 3510966 -y cfrif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_cFRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_plus_coverage.bed` (146027)
<pre>
Cumulative cfrif plot completed!

</pre>
Command completed. Elapsed time: 0:00:36. Running peak memory: 3.534GB.  
  PID: 146027;	Command: Rscript;	Return code: 0;	Memory used: 0.527GB

> `cFRiF`	QC_hg38/K562_RNA-seq_40_cFRiF.pdf	cFRiF	QC_hg38/K562_RNA-seq_40_cFRiF.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_FRiF.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R frif -s K562_RNA-seq_40 -z 3099922541 -n 3510966 -y frif --reads -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_FRiF.pdf --bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_3_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_5_UTR_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Enhancer_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Exon_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Intron_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_plus_coverage.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_Promoter_Flanking_Region_plus_coverage.bed` (146073)
<pre>
Cumulative frif plot completed!

</pre>
Command completed. Elapsed time: 0:00:38. Running peak memory: 3.534GB.  
  PID: 146073;	Command: Rscript;	Return code: 0;	Memory used: 0.381GB

> `FRiF`	QC_hg38/K562_RNA-seq_40_FRiF.pdf	FRiF	QC_hg38/K562_RNA-seq_40_FRiF.png	PEPPRO	_OBJ_

### Calculate mRNA contamination (02-27 10:01:25) elapsed: 76.0 _TIME_

Missing stat 'mRNA_contamination'
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_exons_sort.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_introns_sort.bed`  

> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_exons.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_exons_sort.bed` (146131,146132)
<pre>
</pre>
Command completed. Elapsed time: 0:00:07. Running peak memory: 3.534GB.  
  PID: 146132;	Command: bedtools;	Return code: 0;	Memory used: 0.08GB  
  PID: 146131;	Command: grep;	Return code: 0;	Memory used: 0.005GB


> `grep -wf /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_keep.txt /scratch/jps3dp/DATA/genomes/hg38/refgene_anno/default/hg38_introns.bed | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt | bedtools sort -i stdin -faidx /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_introns_sort.bed` (146140,146141,146142)
<pre>
</pre>
Command completed. Elapsed time: 0:00:08. Running peak memory: 3.534GB.  
  PID: 146141;	Command: bedtools;	Return code: 0;	Memory used: 0.078GB  
  PID: 146140;	Command: grep;	Return code: 0;	Memory used: 0.005GB  
  PID: 146142;	Command: bedtools;	Return code: 0;	Memory used: 0.092GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_coverage.bed`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_coverage.bed`  

> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_exons_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_coverage.bed` (146150)
<pre>
</pre>
Command completed. Elapsed time: 0:00:17. Running peak memory: 3.534GB.  
  PID: 146150;	Command: bedtools;	Return code: 0;	Memory used: 0.022GB


> `bedtools coverage -sorted -counts -s -a /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/hg38_introns_sort.bed -b /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_sort.bam -g /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/chr_order.txt > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_coverage.bed` (146166)
<pre>
</pre>
Command completed. Elapsed time: 0:00:15. Running peak memory: 3.534GB.  
  PID: 146166;	Command: bedtools;	Return code: 0;	Memory used: 0.029GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.243724)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_rpkm.bed` (146180,146181,146182)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 146180;	Command: awk;	Return code: 0;	Memory used: 0.005GB  
  PID: 146182;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 146181;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_rpkm.bed`  

> `awk -v OFS='	' '{chrom[$4] = $1; if($4!=prev4) {chromStart[$4] = $2} strand[$4] = $6; readCount[$4] += $7; exonCount[$4] += 1; geneSizeKB[$4] += (sqrt(($3-$2+0.00000001)^2)/1000); gene[$4] = $4; chromEnd[$4]=$3; prev4=$4} END { for (a in readCount) { print chrom[a], chromStart[a], chromEnd[a], gene[a], (readCount[a]/7.243724)/geneSizeKB[a], strand[a]}}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_coverage.bed | awk '$5>0' | sort -k4 > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_rpkm.bed` (146185,146186,146187)
<pre>
</pre>
Command completed. Elapsed time: 0:00:01. Running peak memory: 3.534GB.  
  PID: 146185;	Command: awk;	Return code: 0;	Memory used: 0.006GB  
  PID: 146187;	Command: sort;	Return code: 0;	Memory used: 0.003GB  
  PID: 146186;	Command: awk;	Return code: 0;	Memory used: 0.001GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed`  

> `join --nocheck-order -a1 -a2 -j4 /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_introns_rpkm.bed /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exons_rpkm.bed | awk -v OFS='	' 'NF==11 {print $7, $8, $9, $1, ($10/$5), $11}' | sort -k1,1 -k2,2n > /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed` (146193,146194,146195)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 146193;	Command: join;	Return code: 0;	Memory used: 0.001GB  
  PID: 146195;	Command: sort;	Return code: 0;	Memory used: 0.004GB  
  PID: 146194;	Command: awk;	Return code: 0;	Memory used: 0.001GB


> `awk '{print $5}' /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'`

> `mRNA_contamination`	3.83	PEPPRO	_RES_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_mRNA_contamination.pdf`  

> `Rscript /scratch/jps3dp/tools/databio/peppro/tools/PEPPRO.R mrna -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed --annotate` (146201)
<pre>
mRNA contamination plot completed!

</pre>
Command completed. Elapsed time: 0:00:06. Running peak memory: 3.534GB.  
  PID: 146201;	Command: Rscript;	Return code: 0;	Memory used: 0.234GB

> `mRNA contamination`	QC_hg38/K562_RNA-seq_40_mRNA_contamination.pdf	mRNA contamination	QC_hg38/K562_RNA-seq_40_mRNA_contamination.png	PEPPRO	_OBJ_
Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed.gz`  

> `pigz -f -p 4 -f /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/QC_hg38/K562_RNA-seq_40_exon_intron_ratios.bed` (146218)
<pre>
</pre>
Command completed. Elapsed time: 0:00:00. Running peak memory: 3.534GB.  
  PID: 146218;	Command: pigz;	Return code: 0;	Memory used: 0.002GB


### Produce bigWig files (02-27 10:02:19) elapsed: 54.0 _TIME_

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam` (146225)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 146225;	Command: samtools;	Return code: 0;	Memory used: 0.01GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (146230)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_plus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_40_plus_cuttrace_6ubo308t'
Processing with 1 cores...
Discarding 129 chunk(s) of reads: ['chrM', 'chr1_KI270708v1_random', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270711v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr9_KI270719v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr14_KI270726v1_random', 'chr17_KI270730v1_random', 'chr22_KI270731v1_random', 'chr22_KI270732v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270583v1', 'chrUn_KI270587v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270589v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270591v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270746v1', 'chrUn_KI270747v1', 'chrUn_KI270749v1', 'chrUn_KI270751v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_GL000218v1', 'chrEBV']
Keeping 66 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr15_KI270727v1_random', 'chr16_KI270728v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chr22_KI270735v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270435v1', 'chrUn_KI270438v1', 'chrUn_KI270580v1', 'chrUn_KI270582v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270741v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270748v1', 'chrUn_KI270750v1', 'chrUn_KI270754v1', 'chrUn_KI270757v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2']
Reduce step (merge files)...
Merging 66 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_exact_body_0-mer.bw'
Merging 66 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_plus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:17:14. Running peak memory: 3.534GB.  
  PID: 146230;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.851GB

Target to produce: `/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_exact_body_0-mer.bw`,`/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_smooth_body_0-mer.bw`  

> `samtools index /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam` (148372)
<pre>
</pre>
Command completed. Elapsed time: 0:00:04. Running peak memory: 3.534GB.  
  PID: 148372;	Command: samtools;	Return code: 0;	Memory used: 0.011GB


> `/scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py -i /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam -c /scratch/jps3dp/DATA/genomes/hg38/fasta/default/hg38.chrom.sizes -o /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_exact_body_0-mer.bw -w /project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_smooth_body_0-mer.bw -p 2 --variable-step --tail-edge` (148376)
<pre>
Cutting parallel chroms in half to accommodate two tracks.
Registering input file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/aligned_hg38/K562_RNA-seq_40_minus.bam'
Temporary files will be stored in: 'tmp_K562_RNA-seq_40_minus_cuttrace_9vaoxar2'
Processing with 1 cores...
stdin is empty of data
Discarding 125 chunk(s) of reads: ['chrM', 'chr1_KI270709v1_random', 'chr1_KI270710v1_random', 'chr1_KI270712v1_random', 'chr2_KI270715v1_random', 'chr2_KI270716v1_random', 'chr5_GL000208v1_random', 'chr9_KI270717v1_random', 'chr9_KI270718v1_random', 'chr14_KI270723v1_random', 'chr14_KI270724v1_random', 'chr14_KI270725v1_random', 'chr16_KI270728v1_random', 'chr22_KI270731v1_random', 'chr22_KI270735v1_random', 'chr22_KI270736v1_random', 'chr22_KI270737v1_random', 'chr22_KI270738v1_random', 'chr22_KI270739v1_random', 'chrY_KI270740v1_random', 'chrUn_KI270302v1', 'chrUn_KI270304v1', 'chrUn_KI270303v1', 'chrUn_KI270305v1', 'chrUn_KI270322v1', 'chrUn_KI270320v1', 'chrUn_KI270310v1', 'chrUn_KI270316v1', 'chrUn_KI270315v1', 'chrUn_KI270312v1', 'chrUn_KI270311v1', 'chrUn_KI270317v1', 'chrUn_KI270412v1', 'chrUn_KI270411v1', 'chrUn_KI270414v1', 'chrUn_KI270419v1', 'chrUn_KI270418v1', 'chrUn_KI270420v1', 'chrUn_KI270424v1', 'chrUn_KI270417v1', 'chrUn_KI270422v1', 'chrUn_KI270423v1', 'chrUn_KI270425v1', 'chrUn_KI270429v1', 'chrUn_KI270465v1', 'chrUn_KI270435v1', 'chrUn_KI270468v1', 'chrUn_KI270510v1', 'chrUn_KI270509v1', 'chrUn_KI270518v1', 'chrUn_KI270508v1', 'chrUn_KI270516v1', 'chrUn_KI270512v1', 'chrUn_KI270519v1', 'chrUn_KI270522v1', 'chrUn_KI270511v1', 'chrUn_KI270515v1', 'chrUn_KI270507v1', 'chrUn_KI270517v1', 'chrUn_KI270529v1', 'chrUn_KI270528v1', 'chrUn_KI270530v1', 'chrUn_KI270539v1', 'chrUn_KI270538v1', 'chrUn_KI270544v1', 'chrUn_KI270548v1', 'chrUn_KI270587v1', 'chrUn_KI270580v1', 'chrUn_KI270581v1', 'chrUn_KI270579v1', 'chrUn_KI270590v1', 'chrUn_KI270584v1', 'chrUn_KI270582v1', 'chrUn_KI270588v1', 'chrUn_KI270593v1', 'chrUn_KI270329v1', 'chrUn_KI270334v1', 'chrUn_KI270333v1', 'chrUn_KI270335v1', 'chrUn_KI270338v1', 'chrUn_KI270340v1', 'chrUn_KI270336v1', 'chrUn_KI270337v1', 'chrUn_KI270363v1', 'chrUn_KI270364v1', 'chrUn_KI270362v1', 'chrUn_KI270366v1', 'chrUn_KI270378v1', 'chrUn_KI270379v1', 'chrUn_KI270389v1', 'chrUn_KI270390v1', 'chrUn_KI270387v1', 'chrUn_KI270395v1', 'chrUn_KI270396v1', 'chrUn_KI270388v1', 'chrUn_KI270394v1', 'chrUn_KI270386v1', 'chrUn_KI270391v1', 'chrUn_KI270383v1', 'chrUn_KI270393v1', 'chrUn_KI270384v1', 'chrUn_KI270392v1', 'chrUn_KI270381v1', 'chrUn_KI270385v1', 'chrUn_KI270382v1', 'chrUn_KI270376v1', 'chrUn_KI270374v1', 'chrUn_KI270372v1', 'chrUn_KI270373v1', 'chrUn_KI270375v1', 'chrUn_KI270371v1', 'chrUn_KI270448v1', 'chrUn_KI270521v1', 'chrUn_KI270741v1', 'chrUn_GL000226v1', 'chrUn_GL000213v1', 'chrUn_KI270747v1', 'chrUn_KI270748v1', 'chrUn_KI270749v1', 'chrUn_KI270752v1', 'chrUn_KI270753v1', 'chrUn_KI270755v1', 'chrUn_KI270756v1', 'chrUn_KI270757v1', 'chrEBV']
Keeping 70 chunk(s) of reads: ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY', 'chr1_KI270706v1_random', 'chr1_KI270707v1_random', 'chr1_KI270708v1_random', 'chr1_KI270711v1_random', 'chr1_KI270713v1_random', 'chr1_KI270714v1_random', 'chr3_GL000221v1_random', 'chr4_GL000008v2_random', 'chr9_KI270719v1_random', 'chr9_KI270720v1_random', 'chr11_KI270721v1_random', 'chr14_GL000009v2_random', 'chr14_GL000225v1_random', 'chr14_KI270722v1_random', 'chr14_GL000194v1_random', 'chr14_KI270726v1_random', 'chr15_KI270727v1_random', 'chr17_GL000205v2_random', 'chr17_KI270729v1_random', 'chr17_KI270730v1_random', 'chr22_KI270732v1_random', 'chr22_KI270733v1_random', 'chr22_KI270734v1_random', 'chrUn_KI270442v1', 'chrUn_KI270466v1', 'chrUn_KI270467v1', 'chrUn_KI270438v1', 'chrUn_KI270583v1', 'chrUn_KI270589v1', 'chrUn_KI270591v1', 'chrUn_KI270330v1', 'chrUn_GL000195v1', 'chrUn_GL000219v1', 'chrUn_GL000220v1', 'chrUn_GL000224v1', 'chrUn_KI270743v1', 'chrUn_KI270744v1', 'chrUn_KI270745v1', 'chrUn_KI270746v1', 'chrUn_KI270750v1', 'chrUn_KI270751v1', 'chrUn_KI270754v1', 'chrUn_GL000214v1', 'chrUn_KI270742v1', 'chrUn_GL000216v2', 'chrUn_GL000218v1']
Reduce step (merge files)...
Merging 70 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_exact_body_0-mer.bw'
Merging 70 files into output file: '/project/shefflab/processed/peppro/paper/dev4/results_pipeline/K562_RNA-seq_40/signal_hg38/K562_RNA-seq_40_minus_smooth_body_0-mer.bw'
</pre>
Command completed. Elapsed time: 0:17:13. Running peak memory: 3.534GB.  
  PID: 148376;	Command: /scratch/jps3dp/tools/databio/peppro/tools/bamSitesToWig.py;	Return code: 0;	Memory used: 1.53GB

Starting cleanup: 52 files; 2 conditional files for cleanup

Cleaning up flagged intermediate files. . .

Cleaning up conditional list. . .

### Pipeline completed. Epilogue
*        Elapsed time (this run):  1:13:24
*  Total elapsed time (all runs):  1:35:16
*         Peak memory (this run):  3.5344 GB
*        Pipeline completed time: 2020-02-27 10:36:56
